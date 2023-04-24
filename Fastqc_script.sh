#!/bin/bash

#megnézendő, hogy fq.gz vagy fq fájlok jöttek-e!

FASTQC=FastQC
RESULTS=FastQC

if [ ! -d $FASTQC ]
then
    mkdir -p $FASTQC
fi

 #fastqc raw_data/*/*.fq.gz -t 3 -o $FastQC

for f in raw_data/*/*_2.fq.gz
do
	dir="$(dirname $f)"
	sample="$(basename $dir)"

	if [ ! -d $RESULTS/QC_raw/$sample ]
	then
		mkdir -p $RESULTS/QC_raw/$sample
	fi

	if [[ $f = *"_2.fq.gz" ]]
	then
		tsamp="$(basename $f _2.fq.gz)"
		first="${dir}/${tsamp}_1.fq.gz"
		second="${dir}/${tsamp}_2.fq.gz"

		if [ -e $first ]
		then
			fastqc $f -t 3 -o $FastQC
			#$FASTQC -o $RESULTS/QC_raw/$sample -t 4 -f fastq $first $second
		fi
	fi
done
