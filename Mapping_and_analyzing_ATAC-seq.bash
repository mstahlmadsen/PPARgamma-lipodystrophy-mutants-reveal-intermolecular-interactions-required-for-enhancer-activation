#!/bin/bash

# Mapping with STAR to mm10, paired end data:
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
	
# Keep primary alignment:
	for SAMFILE in $(ls -d *Aligned.out.sam); do
	echo "Keep only primary alignment for " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 83 || $2 == 99 || $2 == 163 || $2 == 147 { print $0 }' $SAMFILE | cat Header - > ${SAMFILE/.sam/.primary.sam}
	rm Aligned.out.sam
	rm Header 
	done
	
# Determine ducblication % using Picard.
	for SAMFILE in $(ls -d *primary.sam); do
	echo "Sort and deduplicate " $SAMFILE
	# Sort sam-file and determine duplication%
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar SortSam I=$SAMFILE O=$SAMFILE.sorted.sam SO=queryname 2>/dev/null
	java -jar /data/bin/picard/picard-tools-2.5.0/picard.jar MarkDuplicates I=$SAMFILE.sorted.sam O=$SAMFILE.dedup.sam M=$SAMFILE.dedup.metric OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=TRUE 2>/dev/null
	DupPercent=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 9)
	DuplicationPercent=$(echo "scale=2; $DupPercent*100" | bc)
	OpticalDuplicates=$(head $SAMFILE.dedup.metric -n 8 | tail -n 1 | cut -f 8)
	echo -e $SAMFILE"\t"$DuplicationPercent >> Picard_PE_QC.txt
	done
	rm $SAMFILE.sorted.sam
	rm $SAMFILE.dedup.metric
	rm $SAMFILE
	
# Split sam file in two files; NFR (under 120 bp) and MNCL (above 120 bp):
	for VAR in $(ls -d *.sam.dedup.sam); do
	echo "Selecting NFR and MNCL for  " $VAR
	samtools view -h $VAR | awk '$9 <= 120 && $9 >= -120 { print $0 }' > tmp.sam
	samtools view -H $VAR > HEADER
	cat HEADER tmp.sam > ${VAR/.sam.dedup.sam/.dedup_under120.sam}
	samtools view -h $VAR | awk '$9 < -120 || $9 > 120 { print $0 }' > tmp.sam
	samtools view -H $VAR > HEADER
	cat HEADER tmp.sam > ${VAR/.sam.dedup.sam/.dedup_MNCL.sam}
	done
	rm tmp.sam
	rm HEADER

# Generate pseudo-SE SAM-files, make tagdirectories and capture QC:
	for SAMFILE in $(ls -d *.primary.dedup_under120.sam); do
	echo "Making pseudo-SE SAM file on " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 83 { print $1"\t16\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Rv.sam
	awk '$2 == 99 { print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Fw.sam
	cat Header $SAMFILE.dedup_Rv.sam $SAMFILE.dedup_Fw.sam > $SAMFILE.dedup_SE.sam
	done
	for SAMFILE in $(ls -d *.primary.dedup_MNCL.sam); do
	echo "Making pseudo-SE SAM file on " $SAMFILE
	samtools view -H $SAMFILE > Header
	awk '$2 == 83 { print $1"\t16\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Rv.sam
	awk '$2 == 99 { print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6"\t*\t0\t0\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16  }' $SAMFILE > $SAMFILE.dedup_Fw.sam
	cat Header $SAMFILE.dedup_Rv.sam $SAMFILE.dedup_Fw.sam > $SAMFILE.dedup_SE.sam
	done
	
	for SAMFILE in $(ls -d *.dedup_SE.sam); do
	echo "Make Tagdirectory on " $SAMFILE
	makeTagDirectory ${SAMFILE/.sam.dedup_SE.sam/.TD\/} $SAMFILE -format sam -keepAll -genome mm10 -checkGC 2> $SAMFILE.homer
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
	echo -e $SAMFILE"\t"$Positions"\t"$FiltReads"\t"$SameStrand"\t"$DiffStrand"\t"$FLENGTH"\t"$GCContent"\t"$PeakWidth"\t"$Reliable >> TagDir_ATAC_PE_QC.txt
	done
	VAR=$(ls *.TD/tagCountDistribution.txt)
	for i in $VAR
	do
	Out1=$(awk -v name=$i 'NR==2 { print$0 }' $i)
	Out2=$(awk -v name=$i 'NR==3 { print$0 }' $i)
	Out3=$(awk -v name=$i 'NR==4 { print$0 }' $i)
	Out4=$(awk -v name=$i 'NR==5 { print$0 }' $i)
	echo -e $i"\t"$Out1"\t"$Out2"\t"$Out3"\t"$Out4 >> tbp_ATAC.txt
	done
	
# Identify peaks and calculate FRIP:
	for TAGDIRECTORY in $(ls -d *.TD); do
	echo "Doing Homer FindPeaks on " $TAGDIRECTORY
	findPeaks $TAGDIRECTORY/ -size 500 1> ${TAGDIRECTORY/.TD/.peaks} 
	done
		
grep '# total peaks' *.peaks
grep '# Approximate IP efficiency' *.peaks

# Make bigwig-files for visualisation in the UCSC-browser:
	for TD in $(ls -d *.TD); do
	echo "Making UCSC-files based on " $TD
	makeUCSCfile $TD -o ${TD/.TD/.bigwig} -bigWig chrom_size_mm10.txt -fsize 1e20 > trackInfo.txt
	done
	mkdir UCSC_files
	mv *.bigwig UCSC_files
	mv trackInfo.txt UCSC_files
	
# Convert sam-file to bam-file:
	for SAMFILE in $(ls -d *.dedup_under120.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam/.bam}
	done
	for SAMFILE in $(ls -d *.dedup_MNCL.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam/.bam}
	done
	mkdir Bam_files_PE
	mv *.bam Bam_files_PE
	for SAMFILE in $(ls -d *.sam.dedup_SE.sam); do
	echo "Sam to bamming " $SAMFILE
	samtools view -@ 16 -bh $SAMFILE > ${SAMFILE/.sam.dedup_SE.sam/.dedup_SE.bam}
	done
	mkdir Bam_files_SE
	mv *.dedup_SE.bam Bam_files_SE
	rm *.sam
	echo "Script ran to end!"
	exit 0

#____________________________________________________________
# Annotate individal TD to HA-PPARgWT peaks, extended to 500 bp:
	annotatePeaks.pl PPARgWT_peaks_factor_lbg20kb_blacklisted.bed mm10 -noadj -size -250,250 -d *'ATAC'.TD > TagCounts_WTpeaks_ATAC_mm10.txt
# Annotate individal TD to HA-PPARgWT peaks +/- 1500 bp around peak center in bins af 50 bp :	
	annotatePeaks.pl HA_peaks_35.bed mm10 -nogene -noann -size -1500,1500 -hist 50 -ghist -d *'ATAC'.TD > TagCounts_WTpeaksTH35Ext_ATAC_hist_adj_mm10.txt
