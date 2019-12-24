#!/bin/bash

MAXNONSNP=2
DIFF_AS=5
THREADS=1
CHEAP=10

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-a|--annovarfile)
		ANNOVARFILE="$2"
		shift
		shift
		;;
		-b|--bamfile)
		BAMFILE="$2"
		shift
		shift
		;;
		-r|--reference)
		REFERENCE="$2"
		shift
		shift
		;;
		-o|--output-file)
		OUTPUT="$2"
		shift
		shift
		;;
		-s|--snp-database)
		SNPDB="$2"
		shift
		shift
		;;
		-m|--max-non-snp)
		MAXNONSNP="$2"
		shift
		shift
		;;
		-d|--diff-alignment-score)
		DIFF_AS="$2"
		shift
		shift
		;;
		-t|--threads)
		THREADS="$2"
		shift
		shift
		;;
		-c|--current-heapsize)
		CHEAP="$2"
		shift
		shift
		;;
		-h|--help)
		HELP="TRUE"
		shift
		;;
	esac
done
set -- "${POSITIONAL[@]}"

ADD_PARAM=''

if [ ! -z "$HELP" ]
then
	java -jar /code/AdditionalBAMStatistics/additionalBAMStatistics.jar --help
	exit 0
fi

if [ -z "$ANNOVARFILE" ]
then
	echo "Please provide an ANNOVAR-annotated file."
	exit -1
fi

if [ ! -f "$ANNOVARFILE" ]
then
	echo "The ANNOVAR file provided does not exist."
	exit -2
fi

if [ -z "$BAMFILE" ]
then
	echo "Please provide a BAM files."
	exit -1
fi

if [ ! -f "$BAMFILE" ]
then
	echo "The BAM file does not exist."
	exit -2
fi

if [ -z "$REFERENCE" ]
then
	echo "Please provide the reference FASTA file."
	exit -1
fi

if [ ! -f "$REFERENCE" ]
then
	echo "The reference FASTA file does not exist."
	exit -2
fi

if [ ! -z "$OUTPUT" ]
then
	([ -e "$OUTPUT" ] || touch "$OUTPUT") && [ ! -w "$OUTPUT" ] && echo "Cannot write to the provided output file" && exit -1
	ADD_PARAM=$(echo "$ADD_PARAM --output-file $OUTPUT")
fi

if [ ! -z "$SNPDB" ]
then
	if [ ! -f "$SNPDB" ]
	then
		echo "The SNP database file does not exist."
		exit -2
	else
		ADD_PARAM=$(echo "$ADD_PARAM --snp-database $SNPDB")
	fi
fi

regex='^[0-9]+$'

if ! [[ $MAXNONSNP =~ $regex ]]
then
	echo "Please provide an integer value for the maximum number of mismatches not reported in the SNP database per read."
	exit -3
else
	ADD_PARAM=$(echo "$ADD_PARAM --max-non-snp $MAXNONSNP")
fi

if ! [[ $DIFF_AS =~ $regex ]]
then
	echo "Please provide an integer values for the maximum difference between the current and the alternative alignment scores."
	exit -3
else
	ADD_PARAM=$(echo "$ADD_PARAM --difference-alignment-scores $DIFF_AS")

fi

if ! [[ $CHEAP =~ $regex ]]
then
	echo "Please provide an integer value for the maximum JAVA heap memory heap size."
	exit -3
fi 

if ! [[ $THREADS =~ $regex ]]
then
	echo "Please provide an integer value for the number of threads."
	exit -3
else
	ADD_PARAM=$(echo "$ADD_PARAM --threads $THREADS")
fi

java -Xmx${CHEAP}G -jar /lustre/scratch116/casm/cgp/users/ms44/tools/tmp/AdditionalBAMStatistics/additionalBamStatistics.jar --input-annovar-file $ANNOVARFILE --input-bam-file $BAMFILE --reference $REFERENCE $ADD_PARAM 

if [[ $? -ne 0 ]]
then
       	echo -e 'Error: An error has been encountered with AnnotateBamStatistics'
       	exit -4
fi
	
exit 0
