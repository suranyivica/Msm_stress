import sys
import os
import subprocess

#command line arguments
args = list(sys.argv)

f = open(args[1], "r")
f1 = f.readlines()

outfile = open(args[2], "w")


header = f1[0].split('\t')


pos = header[1]
ref = header[2]
mut = header[3]

bam_names = []

for stress in range(4, 73):
    bam_names.append(header[stress])



#print(bam_names)

outfile.write("\t".join(["pos", "ref", "mut"]) + "\t" + "\t".join(bam_names)) 

for i in f1[1:]:
    columns = i.split('\t')
    if float(columns[4]) > 0:
        pass
    elif float(columns[5]) >= 0.06:
        counter = 0
        for stress in range(6, 73):
            try:
                if float(columns[stress]) > 0.0:
                    break
                else:
                    counter += 1
            except ValueError:
                counter += 1
        if counter == 67:     
            outfile.write('\t'.join([columns[1], columns[2], columns[3]]) + "\t" + "\t".join(columns[4:73]) + "\t" + "\n")
            #print(columns[4:73])
	
outfile.close()