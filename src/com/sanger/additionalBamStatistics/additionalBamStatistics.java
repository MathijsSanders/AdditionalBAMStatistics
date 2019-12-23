package com.sanger.additionalBamStatistics;

import java.io.*;
import java.util.*;
import com.beust.jcommander.*;
import com.beust.jcommander.validators.PositiveInteger;

public class additionalBamStatistics {
	private static String versionNumber = "0.3";
	@Parameter
	private List<String> parameters = new ArrayList<String>();
	
	@Parameter(names = "--input-annovar-file", description = "Input ANNOVAR file for further annotating", required = true, converter = FileConverter.class, order=0)
	public File input_annovar_file = null;
	
	@Parameter(names = "--input-bam-file", description = "Input BAM file for processing", required = true, converter = FileConverter.class, order=1)
	public File input_bam_file = null;
	
	@Parameter(names = "--reference", description = "Reference FASTA sequence file used for alignment", required = true, converter = FileConverter.class, order=2)
	public File reference = null;
	
	@Parameter(names = "--output-file", description = "Output file for amended ANNOVAR file (Default: standard out)", order=3)
	public File output_file = null;
	
	@Parameter(names = "--snp-database", description = "SNP database" , converter = FileConverter.class, order=4)
	public File snpdb = null;
	
	@Parameter(names = "--max-non-snp", description = "Maximum number of mismatches in general or not reported in SNP database per read (Dependent on whether SNP database is provided).", validateWith = PositiveInteger.class, order=6)
	public Integer snpMisMatch = 2;
	
	@Parameter(names = "--difference-alignment-scores", description = "Difference threshold between current alignment score and alternative alignment score.", validateWith = PositiveInteger.class, order=7)
	public Integer diffScore = 5;
	
	@Parameter(names = "--threads", description = "Number of threads", validateWith = PositiveInteger.class, order=8)
	public Integer threads = 1;
	
	@Parameter(names = {"--help","-help"}, help = true, description = "Get usage information")
	private boolean help;
	
	@Parameter(names = {"--version","-version"}, description = "Get current version")
	private boolean version;
	public static void main(String[] args) {
		additionalBamStatistics abs  = new additionalBamStatistics();
		JCommander jCommander = new JCommander(abs);
		jCommander.setProgramName("additionalBamStatistics.jar");
		JCommander.newBuilder().addObject(abs).build().parse(args);
		if(abs.version) {
			System.out.println("Calculate BAM statistics from ANNOVAR file: " + versionNumber);
			System.exit(0);
		}
		else if(abs.help) {
			jCommander.usage();
			System.exit(0);
		} else {
			new additionalBamStatisticsCore(abs.input_annovar_file, abs.input_bam_file, abs.reference, abs.snpdb, abs.output_file, abs.diffScore, abs.snpMisMatch, abs.threads);
		}
	}
}