#!/bin/bash

# Mapping with STAR:
	for KEY in $(ls -d *.fastq.gz | awk 'BEGIN{FS = "_"} {print $1}' | sort -u ); do
	FILES=$(ls -d $KEY*.fastq.gz)
	echo "Mapping $KEY, with files: $FILES"
	STAR --genomeLoad LoadAndRemove --genomeDir /data/Genomes/mouse/mm10/star/index101bp/ --runThreadN 10 --readFilesCommand zcat --readFilesIn ${KEY}*.fastq.gz --outFilterMismatchNmax 2 --alignIntronMax 1 --outSAMstrandField none --outSJfilterIntronMaxVsReadN 0 --outFilterMatchNmin 25  --outFileNamePrefix ${KEY}_star.
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
	
# Select primary alignments:
	for SAMFILE in $(ls -d *Aligned.out.sam); do
	echo "Keep only primary alignment for " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 83 || $2 == 99 || $2 == 163 || $2 == 147 { print $0 }' $SAMFILE | cat Header - > ${SAMFILE/.sam/.primary.sam}
	rm Aligned.out.sam
	rm Header 
	done

# Sort sam-file and determine duplication % using Piccard. 
	for SAMFILE in $(ls -d *primary.sam); do
	echo "Sort and deduplicate " $SAMFILE 
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar SortSam I=$SAMFILE O=$SAMFILE.sorted.sam SO=queryname 2>/dev/null
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar MarkDuplicates I=$SAMFILE.sorted.sam O=$SAMFILE.dedup.sam M=$SAMFILE.dedup.metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=TRUE 2>/dev/null
	DupPercent=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 9)
	DuplicationPercent=$(echo "scale=2; $DupPercent*100" | bc)
	OpticalDuplicates=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 8)
	rm $SAMFILE.sorted.sam
	rm $SAMFILE.dedup.metric
# Determine fragment length 
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar SortSam I=$SAMFILE.dedup.sam O=$SAMFILE.dedup.sorted.sam SO=coordinate 2>/dev/null
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar CollectInsertSizeMetrics I=$SAMFILE.dedup.sorted.sam O=$SAMFILE.insert H=$SAMFILE.histo 2>/dev/null
	FLENGTH=$(head $SAMFILE.insert -n 8 | tail -n 1 | cut -f 1) 
	echo -e $SAMFILE"\t"$DuplicationPercent"\t"$FLENGTH >> Picard_PE_QC.txt
	rm $SAMFILE
	done	
	
# Generate pseudo-SE SAM files
	for SAMFILE in $(ls -d *.primary.sam.dedup.sam); do
	echo "Make pseudo-SE Sam files based on " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 83 { print $1"\t16\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Rv.sam
	awk '$2 == 99 { print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Fw.sam
	cat Header $SAMFILE.dedup_Rv.sam $SAMFILE.dedup_Fw.sam > $SAMFILE.dedup_SE.sam
	done
	
# Make tagdirectory and capture output
	for SAMFILE in $(ls -d *.dedup_SE.sam); do
	echo "Make Tagdirectory on " $SAMFILE
	makeTagDirectory ${SAMFILE/.sam.dedup.sam.dedup_SE.sam/.dedup.TD\/} $SAMFILE -format sam -keepAll -fragLength 190 -genome mm10 -checkGC -tbp 1 2> $SAMFILE.homer
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
	echo -e $SAMFILE"\t"$Positions"\t"$FiltReads"\t"$SameStrand"\t"$DiffStrand"\t"$FLENGTH"\t"$GCContent"\t"$PeakWidth"\t"$Reliable >> TagDir_PE_QC.txt
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
	
# Convert sam-files to bam-files:
	for SAMFILE in $(ls -d *.sam.dedup.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam.dedup.sam/.dedup.bam}
	done
	mkdir Bam_files_PE
	mv *.dedup.bam Bam_files_PE
	
	for SAMFILE in $(ls -d *.sam.dedup_SE.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam.dedup_SE.sam/.dedup_SE.bam}
	done
	mkdir Bam_files_SE
	mv *.dedup_SE.bam Bam_files_SE
	done
	echo "Script ran to end!"
	exit 0	
	
# ____________________________________________________________	
# Annotate individual TD to HA-PPARgWT peak-file:
	annotatePeaks.pl PPARgWT_peaks_factor_lbg20kb_blacklisted.bed mm10 -noadj -size -250,250 -d *'MED1'.TD > TagCounts_WTpeaks_MED1_mm10.txt # to feed into R
	annotatePeaks.pl PPARgWT_peaks_factor_lbg20kb_blacklisted.bed mm10 -noann -nogene -size -250,250 -d *'MED1'.TD > TagCounts_WTpeaks_adj_MED1_mm10.txt # to get size factor to feed into R
	annotatePeaks.pl HA_peaks_35.bed  mm10 -nogene -noann -size -1500,1500 -hist 50 -ghist -d *.TD > TagCounts_WTpeaksTH35Ext_MED1_hist_adj_mm10.txt
	