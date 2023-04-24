#Usage: python3 mutation_analysis.py [allele_freq_file] [output file] 

import sys

args = list(sys.argv)

with open(args[1], "r") as f:
    f1 = f.readlines()

#f = open("allele_freq/getaf_output_Comb3_filtered.txt", "r")


    TA_CG = 0
    CG_TA = 0
    TA_AT = 0
    CG_AT = 0
    TA_GC = 0
    CG_GC = 0
    ins = 0
    deletion = 0
    for i in f1:
        columns = i.split('\t')
        if columns[2] == "T" and columns[3] == "C" or columns[2] =="A" and columns[3] == "G":
            if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            TA_CG += adduktum
        elif columns[2] == "C" and columns[3] == "T" or columns[2] =="G" and columns[3] == "A":
            if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            CG_TA += adduktum
        elif columns[2] == "T" and columns[3] == "A" or columns[2] =="A" and columns[3] == "T":
            if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            TA_AT += adduktum
        elif columns[2] == "C" and columns[3] == "A" or columns[2] =="G" and columns[3] == "T":
            if float(columns[5]) >=0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            CG_AT += adduktum
        elif columns[2] == "T" and columns[3] == "G" or columns[2] =="A" and columns[3] == "C":
            if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            TA_GC += adduktum
        elif columns[2] == "C" and columns[3] == "G" or columns[2] =="G" and columns[3] == "C":
            if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                adduktum = 1
            elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                adduktum = 2
            elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                adduktum = 3
            elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                adduktum = 4
            elif float(columns[5]) >= 0.9:
                adduktum = 5
            else:
                adduktum = 0
            CG_GC += adduktum
        else: 
            if len(columns[2]) < len(columns[3]):
                if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                    adduktum = 1
                elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                    adduktum = 2
                elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                    adduktum = 3
                elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                    adduktum = 4
                elif float(columns[5]) >= 0.9:
                    adduktum = 5
                else:
                    adduktum = 0
                ins += adduktum
            elif len(columns[2]) > len(columns[3]):
                if float(columns[5]) >= 0.06 and float(columns[5]) < 0.3:
                    adduktum = 1
                elif float(columns[5]) >= 0.3 and float(columns[5]) < 0.5:
                    adduktum = 2
                elif float(columns[5]) >= 0.5 and float(columns[5]) < 0.7:
                    adduktum = 3
                elif float(columns[5]) >= 0.7 and float(columns[5]) < 0.9:
                    adduktum = 4
                elif float(columns[5]) >= 0.9:
                    adduktum = 5
                else:
                    adduktum = 0
                deletion += adduktum

    output = open(args[2], 'w+')

    output.write('TA to CG' + "\t" + str(TA_CG) + "\n")
    output.write('CG to TA' + "\t" + str(CG_TA) + "\n")
    output.write('TA to AT' + "\t" + str(TA_AT) + "\n")
    output.write('CG to AT' + "\t" + str(CG_AT) + "\n")
    output.write('TA to GC' + "\t" + str(TA_GC) + "\n")
    output.write('CG to GC' + "\t" + str(CG_GC) + "\n")
    output.write('insertion' + "\t" + str(ins) + "\n")
    output.write('deletion' + "\t" + str(deletion) + "\n")
	
	
