#$ -S /bin/bash
#$ -N "unique_hits"
#$ -o /vertessynas/suranyi/2021_02_19_All_hits/unique_hits/unique_hits/output
#$ -e /vertessynas/suranyi/2021_02_19_All_hits/unique_hits/unique_hits/error
#$ -wd /vertessynas/suranyi/2021_02_19_All_hits/unique_hits/unique_hits
#$ -m bea
#$ -M bottgervica@gmail.com
#$ -pe threads 8

DATA=/vertessynas/suranyi/2021_02_19_All_hits/Getaf_for_every_mutations_with_every_strains
Python_file=/vertessynas/suranyi/2021_02_19_All_hits/unique_hits/Filtering_for_unique_mutations.py

export PATH=/bigdisk/programs/bin:$PATH
export PATH=/home/vertessy/suranyi/Python_3.8.5/bin/bin:$PATH

  
for f in $DATA/*all_hits.txt
do        

        sample="$(basename $f)" 
        sample_name="$(cut -d"_" -f1 <<< "$sample")"
        output_unique_hits="${sample_name}_unique_hits.txt"
          

        python3 $Python_file $f $output_unique_hits 
       

done

exit 0
