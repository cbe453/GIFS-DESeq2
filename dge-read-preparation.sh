#!/bin/bash 
# script for running DGE read prep
# this script is meant to be used in conjunction with the group-dge-read-prep.py 
# companion script. Please read the documentation found in the previous script
# before trying to run this one.

set -e

# organizing inputs
GENOTYPE=$1
THREADS=$2
READS_DIR=$3
GENOME=$4
GFF=$5
echo `date`

# check existence of genome and gff files etc.
if [ -e $GENOME.1.ht2 ] 
then
	echo 'Genome file found.'
else
	echo 'Genome file not found. Exiting...'
	exit 1
fi

if [ -e $GFF ] 
then
	echo 'GFF file found.'
else
	echo 'Annotation file (GFF) not found. Exiting...'
	exit 1
fi

if [ -d $READS_DIR ] 
then
	echo 'Reads directory found.'
else
	echo 'Reads directory not found. Exiting...'
	exit 1
fi
	
echo "Sample: $1"
echo "Num threads: $2"
echo "Reads dir: $3"
echo "Genome: $GENOME"
echo "GFF file: $GFF"

# make directory to hold reads from the sample of interest and move
# reads into that directory
#mkdir $GENOTYPE
#mv samples.txt $GENOTYPE
cd $GENOTYPE

if [ -d trim-reads ]
then
	echo "Trim-reads directory already exists."
else
	mkdir trim-reads
fi	

cd trim-reads

for SAMPLE in `cat ../../$GENOTYPE*.txt`
do
	if [ -e $SAMPLE\paired_R1.fastq ]
	then
		echo "Trim files found in trim-reads directory for $SAMPLE\..."
		echo "Skipping $SAMPLE\. Please check that analysis was completed for $SAMPLE"
	else
	    trimmomatic PE -threads $THREADS $READS_DIR/$SAMPLE\R1.fastq $READS_DIR/$SAMPLE\R2.fastq \
			$SAMPLE\paired_R1.fastq $SAMPLE\unpaired_R1.fastq $SAMPLE\paired_R2.fastq \
			$SAMPLE\unpaired_R2.fastq SLIDINGWINDOW:4:25 LEADING:28 TRAILING:28 MINLEN:90

	    if [[ $? -eq 1 ]]; then
		exit 1
	    fi

	    fastqc -t $THREADS $SAMPLE\paired_R1.fastq $SAMPLE\paired_R2.fastq

	    if [[ $? -eq 1 ]]; then
		exit 1
	    fi
	fi
done
cd ../

# directory management and alignment of reads using HiSAT2
if [ -d alignments ]
then
	echo "Alignments directory already exists."
else
	mkdir alignments
fi

cd alignments

for FILE in ../trim-reads/*_paired_R1.fastq
do
	SAMPLE=`basename $FILE | sed -E 's/_paired_R[12].fastq//g'`
	if [ -e $GENOTYPE\_$SAMPLE.sam ] 
	then 
		echo "Alignment file for $SAMPLE is already in alignments directory."
		echo "Skipping $SAMPLE\. Please check that analysis was completed for $SAMPLE"
	else
		echo "Aligning reads to genome for $SAMPLE"
		hisat2 -p $THREADS --dta -q -x $GENOME -1 ../trim-reads/$SAMPLE\_paired_R1.fastq -2 \
		       ../trim-reads/$SAMPLE\_paired_R2.fastq -S $GENOTYPE\_$SAMPLE\.sam

		if [[ $? -eq 1 ]]; then
                    exit 1
		fi
	fi
done

declare -i SAM_THREADS
SAM_THREADS=$THREADS/2
# convert and sort .sam files to .bam
for FILE in ./*.sam
do
	SAMPLE=`basename $FILE | sed -E 's/.sam//g'`
	if [ -e $SAMPLE\.bam ]
	then
		echo "BAM file for $SAMPLE is already in alignments directory."
		echo "Skipping $SAMPLE\. Please check that analysis was completed for $SAMPLE"
	else
		echo "Converting alignments to BAM format"
		samtools view -@ $SAM_THREADS -u $SAMPLE\.sam | samtools sort -@ $SAM_THREADS > $SAMPLE\.bam

		if [[ $? -eq 1 ]]; then
                    exit 1
		fi
	fi
done

cd ../

# directory management for StringTie output
if [ -d stringtie-output ] 
then 
	echo "Stringtie output directory already exists."
else
	mkdir stringtie-output
fi	
	
cd stringtie-output

# Produce gtf files for the later called prepDE.py script
for FILE in ../alignments/*.bam
do
	SAMPLE=`basename $FILE | sed -E 's/.bam//g'`
	if [ -e $SAMPLE\.gtf ]
	then
		echo "GTF file for $SAMPLE is already in stringtie directory."
		echo "Skipping $SAMPLE\. Please check that analysis was completed for $SAMPLE"
	else
		echo "Calling StringTie on $SAMPLE for gene/transcript counts"
		stringtie ../alignments/$SAMPLE\.bam -G $GFF -e -p $THREADS -o $SAMPLE\.gtf

		if [[ $? -eq 1 ]]; then
                    exit 1
		fi
	fi
done

cd ../
echo "Read preparation completed for $GENOTYPE"
echo `date`
