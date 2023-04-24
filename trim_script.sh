#!/bin/bash

DATA=/home/svica/Documents/2019_11_26_Msm_genome_seq
TRIM=/home/svica/Documents/2019_11_26_Msm_genome_seq/Trimmed

TRIMMOMATIC=/home/svica/Downloads/Trimmomatic-0.38/trimmomatic-0.38.jar

if [ ! -d $TRIM ]
then
    mkdir -p $TRIM
fi



for f in $DATA/raw_data/*/*.fq.gz
do
	dir="$(dirname $f)"
	sample="$(basename $dir)"
	first_paired="${sample}_1_paired.fq.gz"
	first_unpaired="${sample}_1_unpaired.fq.gz"
	second_paired="${sample}_2_paired.fq.gz"
	second_unpaired="${sample}_2_unpaired.fq.gz"
	textoutput="${sample}.txt"

	if [[ $f = *"_1.fq.gz" ]]
	then
		echo $f
		tsamp="$(basename $f _1.fq.gz)"
		first="${dir}/${tsamp}_1.fq.gz"
		second="${dir}/${tsamp}_2.fq.gz"

		if [ "$sample" = "Comb3" ] || [ "$sample" = "CycloB" ] || [ "$sample" = "LzdA" ] || [ "$sample" = "LzdB" ] || [ "$sample" = "MockA" ] || [ "$sample" = "MockB" ] || [ "$sample" = "MockC" ]
		then
			java -jar $TRIMMOMATIC PE -phred33 $first $second $TRIM/$first_paired $TRIM/$first_unpaired $TRIM/$second_paired $TRIM/$second_unpaired ILLUMINACLIP:TruSeq2-PE_mod.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee $TRIM/$textoutput
		fi
	fi
done

#TruSeq fájl abban a mappában kell legyen, ahol a script futását elindítom
