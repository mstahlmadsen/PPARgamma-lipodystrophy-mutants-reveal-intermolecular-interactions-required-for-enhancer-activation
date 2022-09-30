#!/bin/bash
# Mapping with STAR and extracting primary reads:	
	for FASTQFILE in $(ls -d *.fq.gz); do
	echo "Mapping" $FASTQFILE
	STAR --genomeLoad LoadAndRemove --genomeDir /data/Genomes/mouse/mm10/star/index101bp/ --readFilesCommand zcat --runThreadN 16 --readFilesIn $FASTQFILE --outFilterMismatchNmax 2 --alignIntronMax 1 --outSAMstrandField none --outSJfilterIntronMaxVsReadN 0 --outFilterMatchNmin 25  --outFileNamePrefix ${FASTQFILE/.fq.gz/_star.}
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
	
	for SAMFILE in $(ls -d *Aligned.out.sam); do
	echo "Keep only primary alignment for " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 0 || $2 == 16 { print $0 }' $SAMFILE | cat Header - > ${SAMFILE/.sam/.primary.sam} 
	rm Aligned.out.sam
	rm Header
	done
	
# Determine duplication:
	for SAMFILE in $(ls -d *primary.sam); do
	echo "Sort and deduplicate " $SAMFILE
	java -jar /data/bin/picard/picard_v2.22.2/picard.jar SortSam I=$SAMFILE O=$SAMFILE.sorted.sam SO=queryname 2>/dev/null
	java -jar /data/bin/picard/picard_v2.22.2/picard.jar MarkDuplicates I=$SAMFILE.sorted.sam O=$SAMFILE.dedup.sam M=$SAMFILE.dedup.metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=TRUE 2>/dev/null
	DupPercent=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 9)
	DuplicationPercent=$(echo "scale=2; $DupPercent*100" | bc)
	OpticalDuplicates=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 8)
	echo -e $SAMFILE"\t"$DuplicationPercent >> Picard_QC.txt
	rm $SAMFILE.sorted.sam
	rm $SAMFILE.dedup.metric
	done

# Make Homer Tag Directories and capture output:
	for SAMFILE in $(ls -d *dedup.sam); do
	echo "Make Tagdirectory on " $SAMFILE
	makeTagDirectory ${SAMFILE/.sam.dedup.sam/.dedup.TD\/} $SAMFILE -format sam -keepAll -genome mm10 -checkGC -tbp 1 2> $SAMFILE.homer
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
	
# Run Homer findPeaks to calculate FRIP:
	for TAGDIRECTORY in $(ls -d *.TD); do
	echo "Doing Homer FindPeaks on " $TAGDIRECTORY
	findPeaks $TAGDIRECTORY/ -size 500 1> ${TAGDIRECTORY/.TD/.TD.peaks} 	
	done
	
# Convert sam-file to bam-file:
	for SAMFILE in $(ls -d *.sam.dedup.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam.dedup.sam/.dedup.bam}
	done
	mkdir Bam_files
	mv *.dedup.bam Bam_files
	rm *.sam
	echo "Script ran to end!"
	exit 0

#________________________________________________________________________________________________________________
## Count ChIP-seq signal within HA-PPARgWT ChIP-seq peaks:
annotatePeaks.pl /data/home/msm/FPLD3/mm10_analysis/HA/HA_peaks_35.bed mm10 -size -250,250 -d *.TD/ > TagCount_PPARgCEBPa_mm10_PPARgWTTH35peaks.txt 
