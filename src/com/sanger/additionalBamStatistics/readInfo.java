package com.sanger.additionalBamStatistics;
import java.util.*;
import htsjdk.samtools.*;
public class readInfo
{
	public boolean varRead = false;
	public int exceedNonSNP = 0;
	public boolean clipped = false;
	public boolean altAlign = false;
	public boolean subExceedAlign = false;
	public boolean smallEdit = false;
	public boolean indel = false;
	public boolean negStrand = false;
	public int variantLoc = -1;
	public boolean isWithinStartFraction = false;
	public boolean containsOtherSomaticVariant = false;
	public HashSet<String> currentSet = null;
	public boolean useForStats = true;
	public boolean hold = false;
	public SAMRecord read = null;
}
