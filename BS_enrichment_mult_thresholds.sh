#!/bin/bash

for REGGENE in $(ls -d *.reg_genes.bed); do
	A=`wc -l < $REGGENE`
	for FACTOR in $(ls -d *.FACT.bed); do
		B=`wc -l < $FACTOR`
		for THRES in 0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 
		do
			slopBed -b $THRES -i $REGGENE -g /data/bin/bedtools/bedtools-v2.29/genomes/mouse.mm10.genome > Gene_region.bed
			intersectBed -wa -a $FACTOR -b Gene_region.bed > Overlap.bed
			C=`wc -l < Overlap.bed`
			echo "scale=4; $C/$A*10000/$B" | bc > ${FACTOR/.FACT.bed/.FACT.$THRES.txt}
			rm Gene_region.bed
			rm Overlap.bed
		done
		echo $FACTOR > FACTOR.name.txt
		paste FACTOR.name.txt *.FACT.*.txt > ${FACTOR/.FACT.bed/.FACT1.$REGGENE.enrich.txt}
		rm *.FACT.*.txt
		rm FACTOR.name.txt
	done
	echo $REGGENE "0 100 10 20 30 40 50 60 70 80 90" > Header.txt
	cat Header.txt *.FACT1.$REGGENE.enrich.txt > Results.$REGGENE.txt
	rm Header.txt
	rm *.enrich.txt
done
cat Results.*.txt > BS.enrichment.results.combined.txt
rm Results.*.txt
echo "I'll bet your binding sites are enriched near regulated genes!"
	