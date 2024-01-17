#$ -S /bin/bash
#$ -N "getaf_for_all_hits"
#$ -o /vertessynas/suranyi/2021_02_19_All_hits/Getaf_for_every_mutations_with_every_strains/output
#$ -e /vertessynas/suranyi/2021_02_19_All_hits/Getaf_for_every_mutations_with_every_strains/error
#$ -wd /vertessynas/suranyi/2021_02_19_All_hits/
#$ -m bea
#$ -pe threads 32

DATA=/vertessynas/suranyi/2020_07_14_unique_check/VCF_copies/rest
BAM=/vertessynas/suranyi/2021_02_19_All_hits/BAM_copies
Python_file=/vertessynas/suranyi/2021_02_19_All_hits/getaf_all_mutations.py

export PATH=/bigdisk/programs/bin:$PATH
export PATH=/home/vertessy/suranyi/Python_3.8.5/bin/bin:$PATH

declare -a BAM_list

for bams in $BAM/*rg.bam
do
        #echo $bams
        bam_name="$(basename $bams)"
        if [[ " ${BAM_list[*]} " == *" $bam_name "* ]] 
        then                 
                : #pass
        else
                BAM_list[${#BAM_list[@]}]=BAM_copies/$bam_name
        fi
        #echo ${BAM_list[@]}    
done
   
for f in $DATA/*_cut
do        

        sample="$(basename $f)" 
        sample_name="$(cut -d"_" -f1 <<< "$sample")"
        output_all_hits="${sample_name}_all_hits.txt"
          

        python3 $Python_file $f $output_all_hits "${BAM_list[@]}" 
       

done
#echo ${BAM_list[@]}

exit 0
