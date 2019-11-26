package com.sanger.additionalBamStatistics;
import java.io.*;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.function.*;
import java.util.concurrent.*;
import java.util.stream.*;

import htsjdk.samtools.*;

public class additionalBamStatisticsCore {
	private Optional<String> header;
	private int exceed = 2;
	private int diffThreshold = 5;
	private final double startFraction = 0.15;
	private final int statsThreshold = 4;
	private final byte qualityThreshold = 25;
	private final int searchWindow = 400;
	public additionalBamStatisticsCore(File inputAnnovar, File inputBam, File referenceFile, File snpDatabase, File outputFile, int diffThreshold, int exceed, int threads) {
		this.diffThreshold = diffThreshold;
		this.exceed = exceed;
		try {
			System.out.println("Loading variants...");
			var variantMap = getVariantList(inputAnnovar);
			System.out.println("Reading reference sequences...");
			getReferenceSequence(variantMap, referenceFile, searchWindow, threads);
			System.out.println("Transforming variants...");
			var variantSet = generateVariantSet(variantMap);
			System.out.println("Loading SNP database...");
			Optional<Map<String, HashSet<String>>> snpdb = Optional.ofNullable((snpDatabase != null) ? readSnpDatabase(snpDatabase) : null);
			System.out.println("Processing variants...");
			for(String chr: variantMap.keySet()) {
				System.out.println("Chromosome " + chr + "...");
				var currentVariants = variantMap.get(chr);
				var currentVariantSet = variantSet.get(chr);
				var currentSNPs = (snpdb.isPresent() ? snpdb.get().get(chr) : null);
				ForkJoinPool forkJoinPool = new ForkJoinPool(threads);
				forkJoinPool.submit(() -> currentVariants.parallelStream().forEach(v -> processVariant(v, currentVariantSet, currentSNPs, inputBam, referenceFile))).get();
			}
			writeResults(variantMap, outputFile);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	private LinkedHashMap<String,ArrayList<variantInfo>> getVariantList(File input) throws IOException {
		Supplier<Stream<String>> streamSupplier = () -> {try{return Files.lines(input.toPath());}catch(IOException e){e.printStackTrace();} return null;};
		header = streamSupplier.get().findFirst();
		return streamSupplier.get().skip(1).map(i -> new variantInfo(i.split("\t"),i)).collect(Collectors.groupingBy(variantInfo::getChromosome,LinkedHashMap::new, Collectors.toCollection(ArrayList::new)));
	}
	private void getReferenceSequence(LinkedHashMap<String, ArrayList<variantInfo>> segregatedVariants, File reference, int searchWindow, int threads) {
		ForkJoinPool forkJoinPool = null;
		try {
			forkJoinPool = new ForkJoinPool(threads);
			forkJoinPool.submit(() -> segregatedVariants.values().parallelStream().forEach(vl -> readSequence(vl, reference, searchWindow))).get();
		} catch (InterruptedException | ExecutionException e) {
			e.printStackTrace();
			System.exit(-3);
		} finally {
			if(forkJoinPool != null)
				forkJoinPool.shutdown();
		}
	}
	private void readSequence(ArrayList<variantInfo> variants, File reference, int searchWindow) {
		var search = variants.stream().map(v -> v.getChromosome() + ":" + (int)(Math.max(1, v.getStart()-searchWindow)) + "-" + (v.getEnd()+searchWindow)).collect(Collectors.toList());
		try {
			String ref = "";
			int track = 0;
			var pb = new ProcessBuilder();
			var command = new ArrayList<String>(Arrays.asList("samtools","faidx",reference.getAbsolutePath()));
			command.addAll(search);
			pb.command(command);
			var proc = pb.start();
			var input = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			String inputLine = null;
			input.readLine();
			while((inputLine = input.readLine()) != null) {
				if(inputLine.startsWith(">")) {
					variants.get(track++).setReference(ref.getBytes());
					ref = "";
				}
				else
					ref += inputLine.toUpperCase();
			}
			variants.get(track).setReference(ref.getBytes());
			input.close();
			var exit = proc.waitFor();
			if(exit != 0) {
				System.out.println("Script has failed");
				System.exit(-4);
			}
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(-4);
		}
	}
	private Map<String, HashSet<String>> generateVariantSet(LinkedHashMap<String, ArrayList<variantInfo>> variantMap) {
		return variantMap.entrySet().stream().map(e -> new AbstractMap.SimpleEntry<String, HashSet<String>>(e.getKey(), e.getValue().stream().map(v -> String.join("_", v.getChromosome(), v.getStart().toString(), v.getRef(), v.getAlt())).collect(Collectors.toCollection(HashSet::new)))).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
	}
	private LinkedHashMap<String,HashSet<String>> readSnpDatabase(File snpDatabase) throws IOException {
		if(isGzipped(snpDatabase))
			return GZIPFiles.lines(snpDatabase.toPath()).filter(s -> !s.startsWith("#")).map(s -> {var b = s.split("\t",6); return Arrays.copyOfRange(b,0,5);}).filter(s -> {if(s[3].length() > 1){return false;}else{var b = s[4].split(","); return Arrays.stream(b).allMatch(com -> com.length() == 1);}}).map(s -> {var b = s[4].split(","); return Arrays.stream(b).map(c -> {var list = Arrays.asList(s[0],s[1],s[3],c);return list.toArray(new String[0]);});}).flatMap(Function.identity()).collect(Collectors.groupingBy(a -> a[0], LinkedHashMap::new, Collectors.mapping(a -> String.join("_",a[0], a[1], a[2], a[3]), Collectors.toCollection(HashSet::new))));
		return Files.lines(snpDatabase.toPath()).filter(s -> !s.startsWith("#")).map(s -> {var b = s.split("\t",6); return Arrays.copyOfRange(b,0,5);}).filter(s -> {if(s[3].length() > 1){return false;}else{var b = s[4].split(","); return Arrays.stream(b).allMatch(com -> com.length() == 1);}}).map(s -> {var b = s[4].split(","); return Arrays.stream(b).map(c -> {var list = Arrays.asList(s[0],s[1],s[3],c);return list.toArray(new String[0]);});}).flatMap(Function.identity()).collect(Collectors.groupingBy(a -> a[0], LinkedHashMap::new, Collectors.mapping(a -> String.join("_",a[0], a[1], a[2], a[3]), Collectors.toCollection(HashSet::new))));
	}
	private boolean isGzipped(File db) {
		int magic = 0;
		try {
			var raf = new RandomAccessFile(db, "r");
			magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
			raf.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-9);
		}
		return magic == GZIPInputStream.GZIP_MAGIC;
	}
	private void processVariant(variantInfo variant, HashSet<String> currentVariantSet, HashSet<String> currentSNPs, File inputBam, File reference) {
		SAMRecord tmpRec = null;
		Cigar cig = null;
		List<CigarElement> el = null;
		ArrayList<readInfo> readList = null;
		var varReadList = new ArrayList<SAMRecord>(200);
		var indelReadList = new ArrayList<SAMRecord>(200);
		var currentSet = new HashSet<String>(100);
		var holdMap = new HashMap<String, readInfo>(100,0.9999f);
		int start = -1;
		var inputSam = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).samRecordFactory(DefaultSAMRecordFactory.getInstance()).open(inputBam);
		SAMRecordIterator it = inputSam.query(variant.getChromosome(), variant.getStart(), variant.getEnd(), false);
		start = ((variant.getStart()-searchWindow) >= 1) ? variant.getStart() - searchWindow : 1;
		var seq = variant.getReference();
		while(it.hasNext()) {
			tmpRec = it.next();
			cig = tmpRec.getCigar();
			el = new ArrayList<CigarElement>(cig.getCigarElements());
			readList = isVarRead(tmpRec, el, seq, start, variant, currentSNPs, startFraction, currentVariantSet, variant.getChromosome(), currentSet, holdMap);
			for(readInfo ri : readList) {
				if(ri.hold)
					continue;
				currentSet = ri.currentSet;
				processRead(ri, variant, varReadList, indelReadList, currentSet, holdMap);
			}
		}
		if(holdMap.keySet().size() > 0)
			for(final readInfo ri : holdMap.values())
				processRead(ri, variant, varReadList, indelReadList, currentSet, holdMap);
		it.close();
		try{inputSam.close();}catch(IOException e) {e.printStackTrace(); System.exit(-5);};
		variant.setVarSomaticReachable(countReachable(varReadList,currentSet));
		variant.setVarReads(varReadList.size());
		variant.setVarReadsUnique(determineUnique(new ArrayList<SAMRecord>(varReadList)));
		variant.setVarReadsUniquePos(determineUniquePos(new ArrayList<SAMRecord>(varReadList)));
		variant.setVarReadsUniqueNeg(determineUniqueNeg(new ArrayList<SAMRecord>(varReadList)));
	}
	private void processRead(readInfo ri, variantInfo dat, ArrayList<SAMRecord> varReadList, ArrayList<SAMRecord> indelReadList, HashSet<String> currentSet, HashMap<String,readInfo> holdMap) {
		currentSet = ri.currentSet;
		if(ri.varRead) {
			varReadList.add(ri.read);
			if(ri.negStrand) {
				dat.setVarNegReads(dat.getVarNegReads() + 1);
				dat.addVarNegLoc(ri.variantLoc);
				dat.addVarNegLocUse(ri.useForStats);
				if(ri.isWithinStartFraction && ri.useForStats)
					dat.setVarNegStartReads(dat.getVarNegStartReads() + 1);
			} else {
				dat.setVarPosReads(dat.getVarPosReads() + 1);
				dat.addVarPosLoc(ri.variantLoc);
				dat.addVarPosLocUse(ri.useForStats);
				if(ri.isWithinStartFraction && ri.useForStats)
					dat.setVarPosStartReads(dat.getVarPosStartReads() + 1);
			}
			if(ri.containsOtherSomaticVariant)
				dat.setVarSomatic(dat.getVarSomatic() + 1);
		}
		if(!ri.varRead && ri.indel)
			indelReadList.add(ri.read);
		if(ri.altAlign) {
			dat.setAltAlignComplete(dat.getAltAlignComplete() + 1);
			if(ri.varRead)
				dat.setAltAlignVar(dat.getAltAlignVar() + 1);
		}
		if(ri.clipped) {
			dat.setClipComplete(dat.getClipComplete() + 1);
			if(ri.varRead)
				dat.setClipVar(dat.getClipVar() + 1);
		}
		if(ri.exceedNonSNP >= exceed) {
			dat.setNMDistanceComplete(dat.getNMDistanceComplete() + 1);
			if(ri.varRead)
				dat.setNMDistanceVar(dat.getNMDistanceVar() + 1);
		}
		if(ri.subExceedAlign) {
			dat.setSubExceedAlignComplete(dat.getSubExceedAlignComplete() + 1);
			if(ri.varRead)
				dat.setSubExceedAlignVar(dat.getSubExceedAlignVar() + 1);
		}
		if(ri.smallEdit)
		{
			dat.setSubSmallEditComplete(dat.getSubSmallEditComplete() + 1);
			if(ri.varRead)
				dat.setSubSmallEditVar(dat.getSubSmallEditVar() + 1);
		}
	}
	private int[] countReachable(ArrayList<SAMRecord> varReadList, HashSet<String> currentSet) {
		if(currentSet.size() == 0)
			return (new int[]{0});
		int[] count = new int[currentSet.size()];
		String[] tokens = null;
		int loc = 0;
		int i;
		for(final SAMRecord read : varReadList) {
			i = 0;
			for(final String current : currentSet) {
				tokens = current.split("_");
				loc = Integer.parseInt(tokens[1]);
				if(read.getAlignmentStart() <= loc && loc <= read.getAlignmentEnd())
					count[i]++;
				i++;
			}
		}
		return count;
	}
	private boolean cigEqualForward(List<CigarElement> cigFirst, List<CigarElement> cigSecond) {
		int i;
		CigarElement telFirst = null;
		CigarElement telSecond = null;
		int track = 0;
		for(i = 0; i < cigFirst.size(); i++) {
			telFirst = cigFirst.get(i);
			if(telFirst.getOperator().toString().equals("S")  || telFirst.getOperator().toString().equals("H"))
				continue;
			if(track == cigSecond.size())
				break;
			telSecond = cigSecond.get(track);
			if(telSecond.getOperator().toString().equals("S")  || telSecond.getOperator().toString().equals("H")) {
				if((track+1) == cigSecond.size())
					break;
				else
					telSecond = cigSecond.get(++track);
			}
			if(!telFirst.getOperator().toString().equals(telSecond.getOperator().toString()) || telFirst.getLength() != telSecond.getLength())
				return false;
			track++;
		}
		return true;
	}
	private boolean cigEqualReverse(List<CigarElement> cigFirst, List<CigarElement> cigSecond) {
		int i;
		CigarElement telFirst = null;
		CigarElement telSecond = null;
		int track = cigSecond.size() - 1;
		for(i = cigFirst.size() - 1; i >= 0; i--) {
			telFirst = cigFirst.get(i);
			if(telFirst.getOperator().toString().equals("S")  || telFirst.getOperator().toString().equals("H"))
				continue;
			if(track == -1)
				break;
			telSecond = cigSecond.get(track);
			if(telSecond.getOperator().toString().equals("S")  || telSecond.getOperator().toString().equals("H")) {
				if((track) == 0)
					break;
				else
					telSecond = cigSecond.get(--track);
			}
			if(!telFirst.getOperator().toString().equals(telSecond.getOperator().toString()) || telFirst.getLength() != telSecond.getLength())
				return false;
			track--;
		}
		return true;
	}
	private int determineUnique(ArrayList<SAMRecord> vList) {
		ArrayList<SAMRecord> unique = new ArrayList<SAMRecord>(50);
		SAMRecord tmp = null;
		SAMRecord current = null;
		Iterator<SAMRecord> it = null;
		while(vList.size() > 0) {
			it = vList.iterator();
			tmp = it.next();
			unique.add(tmp);
			it.remove();
			while(it.hasNext()) {
				current = it.next();
				if(tmp.getReadNegativeStrandFlag()) {
					if(tmp.getAlignmentEnd() == current.getAlignmentEnd() && tmp.getReadNegativeStrandFlag() == current.getReadNegativeStrandFlag())
						if(cigEqualReverse(tmp.getCigar().getCigarElements(), current.getCigar().getCigarElements()))
							it.remove();
				} else {
					if(tmp.getAlignmentStart() == current.getAlignmentStart() && tmp.getReadNegativeStrandFlag() == current.getReadNegativeStrandFlag())
						if(cigEqualForward(tmp.getCigar().getCigarElements(), current.getCigar().getCigarElements()))
							it.remove();
				}
			}
		}
		return unique.size();
	}
	private int determineUniquePos(ArrayList<SAMRecord> vList) {
		ArrayList<SAMRecord> unique = new ArrayList<SAMRecord>(50);
		SAMRecord tmp = null;
		SAMRecord current = null;
		Iterator<SAMRecord> it = null;
		while(vList.size() > 0) {
			it = vList.iterator();
			tmp = it.next();
			if(tmp.getReadNegativeStrandFlag()) {
				it.remove();
				continue;
			} else {
				unique.add(tmp);
				it.remove();
			}
			while(it.hasNext()) {
				current = it.next();
				if(tmp.getAlignmentStart() == current.getAlignmentStart() && tmp.getReadNegativeStrandFlag() == current.getReadNegativeStrandFlag())
					if(cigEqualForward(tmp.getCigar().getCigarElements(), current.getCigar().getCigarElements()))
						it.remove();
			}
		}
		return unique.size();
	}
	private int determineUniqueNeg(ArrayList<SAMRecord> vList) {
		ArrayList<SAMRecord> unique = new ArrayList<SAMRecord>(50);
		SAMRecord tmp = null;
		SAMRecord current = null;
		Iterator<SAMRecord> it = null;
		while(vList.size() > 0) {
			it = vList.iterator();
			tmp = it.next();
			if(!tmp.getReadNegativeStrandFlag()) {
				it.remove();
				continue;
			} else {
				unique.add(tmp);
				it.remove();
			}
			while(it.hasNext()) {
				current = it.next();
				if(tmp.getAlignmentEnd() == current.getAlignmentEnd() && tmp.getReadNegativeStrandFlag() == current.getReadNegativeStrandFlag())
					if(cigEqualReverse(tmp.getCigar().getCigarElements(), current.getCigar().getCigarElements()))
						it.remove();
			}
		}
		return unique.size();
	}
	private void writeResults(LinkedHashMap<String, ArrayList<variantInfo>> variantMap, File outFile) {
		try {
			BufferedWriter output = new BufferedWriter(new FileWriter(outFile));
			StringBuffer buffer = new StringBuffer();
			int counter = 0;
			String reachable = "";
			int[] tmpReach = null;
			int i;
			output.write(header.get() + "\tAlternative_alignment\tAlternative_alignment_var\tSub_Exceed_Align\tSub_Exceed_Align_var\tSub_Small_Edit\tSub_Small_Edit_var\tClipped\tClipped_var\tNonSNP_variants\tNonSNP_variants_var\tVar_pos_reads\tVar_pos_reads_used_stats\tVar_pos_5_prime_15%_reads\tMAD_var_pos_reads\tSD_var_pos_reads\tVar_neg_reads\tVar_neg_reads_used_stats\tVar_neg_5_prime_15%_reads\tMAD_var_neg_reads\tSD_var_neg_reads\tVar_reads\tVar_reads_unique\tVar_reads_unique_pos\tVar_reads_unique_neg\tVar_reads_somatic_variants\tVar_reads_somatic_variants_measurable\n");
			for(final String chr : variantMap.keySet()) {
				var variantList = variantMap.get(chr);
				for(final variantInfo dat: variantList) {
					if(dat.getVarSomaticReachable().length > 1)	{
						tmpReach = dat.getVarSomaticReachable();
						reachable = Integer.toString(tmpReach[0]);
						for(i = 1; i < tmpReach.length; i++)
							reachable += "," + Integer.toString(tmpReach[i]);
					}
					else
						reachable = Integer.toString(dat.getVarSomaticReachable()[0]);	
					buffer.append(dat.getLine() + "\t" + dat.getAltAlignComplete() + "\t" + dat.getAltAlignVar() + "\t" + dat.getSubExceedAlignComplete() + "\t" + dat.getSubExceedAlignVar() + "\t" + dat.getSubSmallEditComplete() + "\t" + dat.getSubSmallEditVar() + "\t" + dat.getClipComplete() + "\t" + dat.getClipVar() + "\t" + dat.getNMDistanceComplete() + "\t" + dat.getNMDistanceVar() + "\t" + dat.getVarPosReads() + "\t" + dat.getCountPosReadsStats() + "\t" + dat.getVarPosStartReads() + "\t" + ((dat.getCountPosReadsStats() >= 2) ? medianAbsoluteDeviation(dat.getVarPosLocs(),dat.getVarPosLocUse()) + "\t" + standardDeviation(dat.getVarPosLocs(),dat.getVarPosLocUse()) : "NA\tNA") + "\t" + dat.getVarNegReads() + "\t" + dat.getCountNegReadsStats() + "\t" + dat.getVarNegStartReads() + "\t" + ((dat.getCountNegReadsStats() >= 2) ? medianAbsoluteDeviation(dat.getVarNegLocs(),dat.getVarNegLocUse()) + "\t" + standardDeviation(dat.getVarNegLocs(),dat.getVarNegLocUse()) : "NA\tNA") + "\t" + dat.getVarReads() + "\t" + dat.getVarReadsUnique() + "\t" + dat.getVarReadsUniquePos() + "\t" + dat.getVarReadsUniqueNeg() + "\t" + dat.getVarSomatic() + "\t" + reachable + "\n");
					counter++;
					if(counter == 100) {
						output.write(buffer.toString());
						output.flush();
						buffer = new StringBuffer();
						counter = 0;
					}
				}
			}
			if(counter > 0) {
				output.write(buffer.toString());
				output.flush();
			}
			output.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	private double standardDeviation(ArrayList<Double> locs, ArrayList<Boolean> useList) {
		double sum = 0.0;
		double tot = 0.0;
		int i;
		for(i = 0; i < useList.size(); i++) {
			if(useList.get(i)) {
				tot += 1.0;
				sum += locs.get(i);
			}
		}
		double mean = sum / tot;
		double sumDev = 0.0;
		for(i = 0; i < useList.size(); i++)
			if(useList.get(i)) 
				sumDev += Math.pow(locs.get(i) - mean, 2.0);
		sumDev *= 1.0/(tot - 1.0);
		return Math.sqrt(sumDev);
	}
	private double medianAbsoluteDeviation(ArrayList<Double> allLocs, ArrayList<Boolean> useList) {
		ArrayList<Double> locs = new ArrayList<Double>(allLocs.size());
		int i;
		for(i = 0; i < useList.size(); i++)
			if(useList.get(i))
				locs.add(allLocs.get(i));
		Double[] val = new Double[locs.size()];
		locs.toArray(val);
		Arrays.sort(val);
		double med = 0.0;
		if(val.length % 2 == 0)
			med = (val[(val.length/2)-1] + val[val.length/2])/2.0;
		else
			med = val[val.length/2];
		double[] dev = new double[val.length];
		int pos = 0;
		for(double loc : val)
			dev[pos++] = Math.abs(loc - med);
		Arrays.sort(dev);
		if(dev.length %2 == 0)
			return ((dev[(dev.length/2)-1] + dev[dev.length/2])/2.0);
		return dev[dev.length/2];
	}
	private ArrayList<readInfo> isVarRead(SAMRecord tmpRec, List<CigarElement> el, byte[] ref, int start, variantInfo dat, HashSet<String> snpMap, double startFraction, HashSet<String> variantSet, String chromosome, HashSet<String> currentSet, HashMap<String, readInfo> holdMap) {
		ArrayList<readInfo> readList = new ArrayList<readInfo>(2);
		readInfo tmp = new readInfo();
		readInfo mate = holdMap.get(tmpRec.getReadName());
		tmp.currentSet = currentSet;
		int trackerRef = tmpRec.getAlignmentStart() - start;
		int trackerSeq = 0;
		double trackerStart = 0.0;
		int vLoc = 0;
		byte[] seqRead = tmpRec.getReadBases();
		byte[] qualityRead = tmpRec.getBaseQualities();
		Iterator<CigarElement> it = el.iterator();
		CigarElement tel = null;
		it = el.iterator();
		int i;
		String search = null;
		int searchLoc = -1;
		String xaInfo = null;
		try {
			xaInfo = tmpRec.getStringAttribute("XA");
			if(xaInfo != null && !xaInfo.equals(""))
				tmp.altAlign = true;
			if(tmpRec.getIntegerAttribute("XS") >= tmpRec.getIntegerAttribute("AS"))
				tmp.subExceedAlign = true;
			if(Math.abs(tmpRec.getIntegerAttribute("AS")-tmpRec.getIntegerAttribute("XS")) <= diffThreshold)
				tmp.smallEdit = true;
			tmp.negStrand = tmpRec.getReadNegativeStrandFlag();
			while(it.hasNext()) {
				tel = it.next();
				if(tel.getOperator().toString().equals("H"))
					tmp.clipped = true;
				else if(tel.getOperator().toString().equals("S")) {
					trackerStart += tel.getLength();
					trackerSeq += tel.getLength();
					tmp.clipped = true;
				}
				else if(tel.getOperator().toString().equals("M")) {
					for(i = 0; i < tel.getLength(); i++) {
						if(seqRead[trackerSeq] == 'N') {
							trackerRef++;
							trackerSeq++;
							trackerStart++;
							vLoc++;
							continue;
						}
						if(ref[trackerRef] != seqRead[trackerSeq]) {
							searchLoc = (start + trackerRef);
							if(searchLoc == dat.getStart() && seqRead[trackerSeq] == dat.getAlt().getBytes()[0]) {
								tmp.varRead = true;
								tmp.variantLoc = vLoc;
								tmp.read = tmpRec;
								tmp.useForStats = (qualityRead[trackerSeq] >= qualityThreshold);
								if(tmpRec.getReadNegativeStrandFlag()) 
									tmp.useForStats = (tmp.useForStats && tmpRec.getAlignmentStart() >= tmpRec.getMateAlignmentStart());
								else
									tmp.useForStats = (tmp.useForStats && tmpRec.getAlignmentStart() <= tmpRec.getMateAlignmentStart());
							} else {
								search = searchLoc + "_" + searchLoc + "_" + (char)ref[trackerRef] + "_" + (char)seqRead[trackerSeq];
								if(snpMap == null || !snpMap.contains(search))
									tmp.exceedNonSNP++;
								search = chromosome + "_" + searchLoc + "_" + (char)ref[trackerRef] + "_" + (char)seqRead[trackerSeq];
								if(variantSet.contains(search)) {
									tmp.containsOtherSomaticVariant = true;
									if(!tmp.currentSet.contains(search))
										tmp.currentSet.add(search);
								}
							}
						}
						trackerRef++;
						trackerSeq++;
						trackerStart++;
						vLoc++;
					}
				}
				else if(tel.getOperator().toString().equals("I")) {
					trackerSeq += tel.getLength();
					vLoc += tel.getLength();
					trackerStart += tel.getLength();
					tmp.indel = true;
				}
				else if(tel.getOperator().toString().equals("D")) {
					trackerRef += tel.getLength();
					tmp.indel = true;
				}
				else {
					System.out.println("You are missing the following operator: " + tel.getOperator().toString());
					System.exit(0);
				}
			}
			if(tmpRec.getReadNegativeStrandFlag())
				tmp.variantLoc = vLoc - tmp.variantLoc;
			tmp.isWithinStartFraction = ((((double)tmp.variantLoc)/trackerStart) <= startFraction);
			tmp.useForStats = (tmp.useForStats && tmp.exceedNonSNP <= statsThreshold);
			if(mate == null && tmp.varRead && tmp.useForStats && !tmpRec.getMateUnmappedFlag() && tmpRec.getAlignmentStart() <= tmpRec.getMateAlignmentStart() && tmpRec.getMateAlignmentStart() <= dat.getStart() && tmpRec.getReferenceIndex() == tmpRec.getMateReferenceIndex()) {
				tmp.hold = true;
				holdMap.put(tmpRec.getReadName(), tmp);
			}
			if(mate != null) {
				mate.hold = false;
				if(mate.varRead && tmp.varRead) {
					if(mate.variantLoc <= tmp.variantLoc) {
						mate.useForStats = true;
						tmp.useForStats = false;
					} else {
						mate.useForStats = false;
						tmp.useForStats = true;
					}
				}
				readList.add(mate);
				holdMap.remove(tmpRec.getReadName());
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		readList.add(tmp);
		return readList;
	}
	
}
