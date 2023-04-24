GENOME=/home/svica/Documents/Msm_reference/Mycobacterium_smegmatis_str_mc2_155.ASM1500v1.dna.toplevel.fa
REF_VCF=/home/svica/Documents/2020_05_22_Msm_genome_seq/VCF/APAW_raw.vcf

VCFDIR=/home/svica/Documents/2020_05_22_Msm_genome_seq/VCF

GATK=~/build/gatk/build/libs/gatk.jar


for f in $VCFDIR/*.vcf
do
	
	sample="$(basename $f _raw.vcf )"
	vcf_shared_output="${VCFDIR}/${sample}_shared.vcf"
	vcf_reduced_output="${VCFDIR}/${sample}_reduced.vcf"
	cut_output="${VCFDIR}/${sample}_cut"

	java -jar $GATK SelectVariants -V $VCFDIR/${sample}_raw.vcf -R $GENOME -conc $REF_VCF -O $vcf_shared_output
	java -jar $GATK SelectVariants -V $VCFDIR/${sample}_raw.vcf -R $GENOME -disc $VCFDIR/${sample}_shared.vcf -O $vcf_reduced_output
	grep -v "#" $vcf_reduced_output |cut -f1-6 > $cut_output

done



