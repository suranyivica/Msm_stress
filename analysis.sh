DATA=/home/svica/Documents/2020_05_22_Msm_genome_seq
GENOME=/home/svica/Documents/Msm_reference/Mycobacterium_smegmatis_str_mc2_155.ASM1500v1.dna.toplevel.fa

SAMDIR=/home/svica/Documents/2020_05_22_Msm_genome_seq/SAM
BAMDIR=/home/svica/Documents/2020_05_22_Msm_genome_seq/BAM
VCFDIR=/home/svica/Documents/2020_05_22_Msm_genome_seq/VCF

PICARD=~/build/picard/build/libs/picard.jar
GATK=~/build/gatk/build/libs/gatk.jar

if [ ! -d $SAMDIR ]
then
    mkdir -p $SAMDIR
fi

if [ ! -d $BAMDIR ]
then
    mkdir -p $BAMDIR
fi

if [ ! -d $VCFDIR ]
then
    mkdir -p $VCFDIR
fi




for f in $DATA/raw_data/*/*_1.fq.gz
do
	dir="$(dirname $f)"
	sample="$(basename $dir)"
	samoutput="${SAMDIR}/${sample}.sam"
	bamoutput="${BAMDIR}/${sample}.bam"
	rgbamoutput="${BAMDIR}/${sample}.rg.bam"
	vcfoutput="${VCFDIR}/${sample}_raw.vcf"

	if [[ $f = *"_1.fq.gz" ]]
	then
		tsamp="$(basename $f _1.fq.gz)"
		first="${dir}/${tsamp}_1.fq.gz"
		second="${dir}/${tsamp}_2.fq.gz"

		if [ -e $first ]
		then
			bowtie2 -p 4 -x $GENOME -1 $first -2 $second -S $samoutput
                       	samblaster -r -i $samoutput | samtools view -buh | samtools sort - > $bamoutput
			java -jar $PICARD AddOrReplaceReadGroups I=$bamoutput O=$rgbamoutput RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
			samtools index $rgbamoutput
			if [ "$sample" = "BB" ]
			then
				java -jar $GATK HaplotypeCaller -R $GENOME -I $rgbamoutput -O $vcfoutput -ploidy 6
			elif [ "$sample" = "C3" ] || [ "$sample" = "C5" ]
			then
				java -jar $GATK HaplotypeCaller -R $GENOME -I $rgbamoutput -O $vcfoutput -ploidy 1
			else
				java -jar $GATK HaplotypeCaller -R $GENOME -I $rgbamoutput -O $vcfoutput -ploidy 5
			fi
		fi
	fi

done

#ezután jön a vcf_work



