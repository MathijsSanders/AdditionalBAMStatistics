package com.sanger.additionalBamStatistics;
import java.util.ArrayList;

public class variantInfo 
{
	private String chromosome = "";
	private Integer start = -1;
	private Integer end = -1;
	private byte[] reference = null;
	private String ref = "";
	private String alt = "";
	private String line = "";
	private int altAlignComplete = 0;
	private int altAlignVar = 0;
	private int nmDistanceComplete = 0;
	private int nmDistanceVar = 0;
	private int clipComplete = 0;
	private int clipVar = 0;
	private int varReads = 0;
	private int varReadsUnique = 0;
	private int varReadsUniquePos = 0;
	private int varReadsUniqueNeg = 0;
	private int subExceedAlignComplete = 0;
	private int subExceedAlignVar = 0;
	private int subSmallEditComplete = 0;
	private int subSmallEditVar = 0;
	private int varPosReads = 0;
	private int varNegReads = 0;
	private int varPosStartReads = 0;
	private int varNegStartReads= 0;
	private ArrayList<Double> posVarLoc = new ArrayList<Double>(100);
	private ArrayList<Boolean> posVarLocUse = new ArrayList<Boolean>(100);
	private ArrayList<Double> negVarLoc = new ArrayList<Double>(100);
	private ArrayList<Boolean> negVarLocUse = new ArrayList<Boolean>(100);
	private int varSomatic = 0;
	private int[] varSomaticReachable = null;
	public variantInfo(String[] info, String line) throws NumberFormatException {
		this.chromosome = info[0];
		this.start = Integer.parseInt(info[1]);
		this.end = Integer.parseInt(info[2]);
		this.ref = info[3];
		this.alt = info[4];
		this.line = line;
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public Integer getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public Integer getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public String getRef() {
		return ref;
	}
	public void setRef(String ref) {
		this.ref = ref;
	}
	public String getAlt() {
		return alt;
	}
	public void setAlt(String alt) {
		this.alt = alt;
	}
	public String getLine() {
		return line;
	}
	public void setAltAlignComplete(int aac) {
		altAlignComplete = aac;
	}
	public void setAltAlignVar(int aav) {
		altAlignVar = aav;
	}
	public int getAltAlignComplete() {
		return altAlignComplete;
	}
	public int getAltAlignVar() {
		return altAlignVar;
	}
	public int getNMDistanceComplete() {
		return nmDistanceComplete;
	}
	public void setNMDistanceComplete(int ndc) {
		nmDistanceComplete = ndc;
	}
	public int getNMDistanceVar() {
		return nmDistanceVar;
	}
	public void setNMDistanceVar(int ndv) {
		nmDistanceVar = ndv;
	}
	public int getClipComplete() {
		return clipComplete;
	}
	public void setClipComplete(int cc) {
		clipComplete = cc;
	}
	public int getClipVar() {
		return clipVar;
	}
	public void setClipVar(int cv) {
		clipVar = cv;
	}
	public int getVarReads() {
		return varReads;
	}
	public void setVarReads(int vr) {
		varReads = vr;
	}
	public int getVarReadsUnique() {
		return varReadsUnique;
	}
	public void setVarReadsUnique(int vru) {
		varReadsUnique = vru;
	}
	public int getVarReadsUniquePos() {
		return varReadsUniquePos;
	}
	public void setVarReadsUniquePos(int vrup) {
		varReadsUniquePos = vrup;
	}
	public int getVarReadsUniqueNeg() {
		return varReadsUniqueNeg;
	}
	public void setVarReadsUniqueNeg(int vrun) {
		varReadsUniqueNeg = vrun;
	}
	public int getSubExceedAlignComplete() {
		return subExceedAlignComplete;
	}
	public void setSubExceedAlignComplete(int seac) {
		subExceedAlignComplete = seac;
	}
	public int getSubExceedAlignVar() {
		return subExceedAlignVar;
	}
	public void setSubExceedAlignVar(int seav) {
		subExceedAlignVar = seav;
	}
	public int getSubSmallEditComplete() {
		return subSmallEditComplete;
	}
	public void setSubSmallEditComplete(int ssec) {
		subSmallEditComplete = ssec;
	}
	public int getSubSmallEditVar() {
		return subSmallEditVar;
	}
	public void setSubSmallEditVar(int ssev) {
		subSmallEditVar = ssev;
	}
	public int getVarPosReads() {
		return varPosReads;
	}
	public void setVarPosReads(int vpr) {
		varPosReads = vpr;
	}
	public int getVarNegReads() {
		return varNegReads;
	}
	public void setVarNegReads(int vnr) {
		varNegReads = vnr;
	}
	public ArrayList<Double> getVarPosLocs() {
		return posVarLoc;
	}
	public void addVarPosLoc(double loc) {
		posVarLoc.add(loc);
	}
	public ArrayList<Boolean> getVarPosLocUse() {
		return posVarLocUse;
	}
	public long getCountPosReadsStats() {
		return posVarLocUse.stream().filter(i -> i == true).count();
	}
	public void addVarPosLocUse(Boolean use) {
		posVarLocUse.add(use);
	}
	public ArrayList<Double> getVarNegLocs() {
		return negVarLoc;
	}
	public long getCountNegReadsStats() {
		return negVarLocUse.stream().filter(i -> i == true).count();
	}
	public void addVarNegLoc(double loc) {
		negVarLoc.add(loc);
	}
	public ArrayList<Boolean> getVarNegLocUse() {
		return negVarLocUse;
	}
	public void addVarNegLocUse(Boolean use) {
		negVarLocUse.add(use);
	}
	public int getVarPosStartReads() {
		return varPosStartReads;
	}
	public void setVarPosStartReads(int vpsr) {
		varPosStartReads = vpsr;
	}
	public int getVarNegStartReads() {
		return varNegStartReads;
	}
	public void setVarNegStartReads(int vnsr) {
		varNegStartReads = vnsr;
	}
	public int getVarSomatic() {
		return varSomatic;
	}
	public void setVarSomatic(int vs) {
		varSomatic = vs;
	}
	public int[] getVarSomaticReachable() {
		return varSomaticReachable;
	}
	public void setVarSomaticReachable(int[] vsr) {
		varSomaticReachable = vsr;
	}
	public byte[] getReference() {
		return reference;
	}
	public void setReference(byte[] reference) {
		this.reference = reference;
	}
}
