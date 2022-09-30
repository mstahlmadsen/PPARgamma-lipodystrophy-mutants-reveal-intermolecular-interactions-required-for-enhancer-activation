#!/bin/bash

# Mapping with STAR
	for FASTQFILE in $(ls -d *.fastq.gz); do
	echo "Mapping" $FASTQFILE
	STAR --genomeLoad LoadAndRemove --genomeDir /data/Genomes/mouse/mm10/star/index101bp/ --readFilesCommand zcat --runThreadN 16 --readFilesIn $FASTQFILE --outFileNamePrefix ${FASTQFILE/_S*_L001_R1_001.fastq.gz/.star_}
	done
	VAR=$(ls *Log.final.out)
		for i in $VAR
		do
		Out1=$(awk -v name=$i 'NR==6 { print$0 }' $i)
		Out2=$(awk -v name=$i 'NR==9 { print$0 }' $i)
		Out3=$(awk -v name=$i 'NR==10 { print$0 }' $i)
		Out4=$(awk -v name=$i 'NR==30 { print$0 }' $i)
	echo -e $i"\t"$Out1"\t"$Out2"\t"$Out3"\t"$Out4 >> Aligned.txt
	done 
	echo "Cleaning up"
	mkdir Log.final.out_files
	mv *Log.final.out Log.final.out_files
	mkdir Fastq_files
	mv *.fastq.gz Fastq_files 
	
# Keep primary alignment:
	for SAMFILE in $(ls -d *Aligned.out.sam); do
	echo "Keep only primary alignment for " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 0 || $2 == 16 { print $0 }' $SAMFILE | cat Header - > ${SAMFILE/.sam/.primary.sam}
	rm Aligned.out.sam
	rm Header 
	done
	
# Determine duplication with Picard:
	for SAMFILE in $(ls -d *primary.sam); do
	echo "Sort and deduplicate " $SAMFILE
	# Sort sam-file and determine duplication%
	java -jar /data/bin/picard/picard_v2.22.2/picard.jar SortSam I=$SAMFILE O=$SAMFILE.sorted.sam SO=queryname 2>/dev/null
	java -jar /data/bin/picard/picard_v2.22.2/picard.jar MarkDuplicates I=$SAMFILE.sorted.sam O=$SAMFILE.dedup.sam M=$SAMFILE.dedup.metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2>/dev/null
	DupPercent=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 9)
	DuplicationPercent=$(echo "scale=2; $DupPercent*100" | bc)
	OpticalDuplicates=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 8)
	rm $SAMFILE.dedup.metric
	echo -e $SAMFILE"\t"$DuplicationPercent >> Picard_SE_QC.txt
	done
	
# Make Homer Tag Directories and extract QC from TD:
	for SAMFILE in $(ls -d *sorted.sam); do
	echo "Make Tagdirectory on " $SAMFILE
	makeTagDirectory ${SAMFILE/.sam.sorted.sam/.TD\/} $SAMFILE -format sam -keepAll -genome mm10 -checkGC 2> $SAMFILE.homer
	Positions=$(head $SAMFILE.homer -n 11 | tail -n 1 | awk '{ print $4 }')
	FiltReads=$(head $SAMFILE.homer -n 10 | tail -n 1 | awk '{ print $4 }')
	SameStrand=$(cat $SAMFILE.homer | grep 'Same strand' | awk '{ print $5 }')
	DiffStrand=$(cat $SAMFILE.homer | grep 'Diff strand' | awk '{ print $5 }')
	FLENGTH=$(cat $SAMFILE.homer | grep 'Current Fragment length' | awk '{ print $5 }')
	GCContent=$(cat $SAMFILE.homer | grep 'Avg Fragment GC%' | awk '{ print $5 }')
	PeakWidth=$(cat $SAMFILE.homer | grep 'Peak Width Estimate' | awk '{ print $4 }')
	Reliable=$(cat $SAMFILE.homer | grep 'No reliable estimate' | wc -l)
	if [ "$Reliable" = "1" ]
	then
	PeakWidth="Unreliable"
	fi
	echo -e $SAMFILE"\t"$Positions"\t"$FiltReads"\t"$SameStrand"\t"$DiffStrand"\t"$FLENGTH"\t"$GCContent"\t"$PeakWidth"\t"$Reliable >> TagDir_SE_QC.txt
	done
	VAR=$(ls *.TD/tagCountDistribution.txt)
	for i in $VAR
	do
	Out1=$(awk -v name=$i 'NR==2 { print$0 }' $i)
	Out2=$(awk -v name=$i 'NR==3 { print$0 }' $i)
	Out3=$(awk -v name=$i 'NR==4 { print$0 }' $i)
	Out4=$(awk -v name=$i 'NR==5 { print$0 }' $i)
	echo -e $i"\t"$Out1"\t"$Out2"\t"$Out3"\t"$Out4 >> tbp.txt
	done

# Make bigwig-files for visualisation in the UCSC-browser:
	for TD in $(ls -d *.TD); do
	echo "Making UCSC-files based on " $TD
	makeUCSCfile $TD -o ${TD/.TD/.bigwig} -bigWig chrom_size_mm10.txt -fsize 1e20 > trackInfo.txt
	done
	mkdir UCSC_files
	mv *.bigwig UCSC_files
	mv trackInfo.txt UCSC_files

# Convert sam to bam-file:
	for SAMFILE in $(ls -d *.sorted.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam.sorted.sam/.sorted.bam}
	done
	mkdir Bam_files
	mv *.sorted.bam Bam_files
	rm *.sam
	
# annotate reads to exons
	analyzeRepeats.pl rna mm10 -count exons -condenseGenes -noCondensing -noadj -d *.TD > Counts_RNA_FPLD3_mm10.txt
	exit 0