import sys
import os
import subprocess

try:
    import pysam
except ImportError:
    import pip
    pip.main(["install", 'pysam'])

# command line arguments
args = list(sys.argv)
bams = args[5:]
if len(args) == 0:
    sys.exit("Usage: python3 getaf.py [position file] [output file for unique mutations] [output file for low frequency unique mutations] [output file for not unique mutations] [bam files]")
if not all(x.endswith("bam") for x in bams):
    sys.exit("Usage: python3 getaf.py [position list] [output file for unique mutations] [output file for low frequency unique mutations] [output file for not unique mutations] [bam files]" + " ".join([bams]))

j = int(subprocess.check_output(['wc', '-l', args[1]]).split()[0])


bam_h = list(map(os.path.basename, bams))
bam_names = []
for elements in bam_h:
    names_parts = elements.split(".")
    bam_names.append(names_parts[0])
    
outfile_unique = open(args[2], "w")
outfile_unique.write("\t".join(["chromosome", "pos", "ref", "mut", "score"]) + "\t" + "\t".join(bam_names) + "\n")

outfile_low_freq_unique = open(args[3], "w")
outfile_low_freq_unique.write("\t".join(["chromosome", "pos", "ref", "mut", "score"]) + "\t" + "\n")


outfile_not_unique_mutations = open(args[4], "w")
outfile_not_unique_mutations.write("\t".join(["chromosome", "pos", "ref", "mut", "score"]) + "\t" + "\t".join(bam_names) + "\n")

linecounter = 0
with open(args[1], "r") as f:
    for line in f:
        linecounter = linecounter + 1
        words = line.strip().split("\t")
        chrom = words[0]
        pos   = words[1]
        ref   = words[3]
        mut   = words[4]
        score = words[5]
        if chrom not in pysam.AlignmentFile(bams[0], "rb").references:
            sys.exit("Error at line" +
            str(linecounter) +
            ": some chromosome names are not in the reference (eg. chr1 instead of 1)")
        values=[]
        coverages=[]
        # what kind of mutation: SNV, del, ins
        if len(mut) == len(ref):
            muttype = "snv"
        elif len(mut) < len(ref):
            muttype = "del"
        else:
            muttype = "ins"
        for b in bams:
            bam = pysam.AlignmentFile(b, "rb")
            counter = 0
            total = 0
            for pileupcolumn in bam.pileup(str(chrom), int(pos) - 1, int(pos), nostepper = "samtools", min_base_quality = 0):
                
                if pileupcolumn.pos == int(pos) - 1:
                    total = pileupcolumn.nsegments
                    #print('position ' + pos)
                    #print('total' + str(total))
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.query_position is not None:
                            # if the mutation type is SNV, I will just extract the actual sequence at the position, and check if it is the same as mut
                            # if it is indel, I will use the CIGAR string
                            if muttype == "snv":
                                if pileupread.alignment.query_sequence[pileupread.query_position] == mut:
                                    counter = counter + 1
                            if muttype == "del":
                                #print('I am in the del side')
                                p = pileupread.alignment.reference_start
                                #if pileupread.alignment.query_name == 'A00808:22:HJ7LMDSXX:3:1371:5421:16423':
                                #print('p' + str(p))
                                #print(pileupread.alignment.query_name)
                                if p == int(pos):
                                    counter = counter + 1
                                for ce in pileupread.alignment.cigartuples:
                                    #print(ce)
                                    
                                    if ce[0] in [0, 2]:
                                        p = p + ce[1]
                                        #print('p' + str(p))
                                        if p == int(pos):
                                            #print('I am in the pos side')
                                            counter = counter + 1
                                            #print(counter)
                                
                            if muttype == "ins":
                                p = pileupread.alignment.reference_start
                                for ce in pileupread.alignment.cigartuples:
                                    if p == int(pos) and ce[0] == 1:
                                        counter = counter + 1
                                    if ce[0] in [0, 2]:
                                        p = p + ce[1]
                                if p == pos:
                                    counter = counter + 1
                    #print('counter ' + str(counter))
            if total == 0:
                values.append("NA")
            else:
                values.append(str(round(counter / total, 3)))
                coverages.append(str(total))
        common_counter = 0
        
            
            
        if float(values[0] >0 or float(values[2]) > 0 or float(values[3]) > 0 or float(values[4]) > 0 or float(values[5]) > 0 or float(values[6]) > 0 or float(values[7]) > 0 or float(values[8]) > 0 or float(values[9]) > 0 or float(values[10]) > 0 or float(values[11]) > 0 or float(values[12]) > 0 or float(values[13]) > 0 or float(values[14]) > 0 or float(values[15]) > 0 or float(values[16]) > 0 or float(values[17]) > 0 or float(values[18]) > 0 or float(values[19]) > 0 or float(values[20]) > 0 or float(values[21]) > 0 or float(values[22]) > 0 or float(values[23]) > 0 or float(values[24]) > 0 or float(values[25]) > 0 or float(values[26]) > 0 or float(values[27]) > 0 or float(values[28]) > 0 or float(values[29]) > 0 or float(values[30]) > 0 or float(values[31]) > 0 or float(values[32]) > 0 or float(values[33]) > 0 or float(values[34]) > 0 or float(values[35]) > 0 or float(values[36]) > 0 or float(values[37]) > 0 or float(values[38]) > 0 or float(values[39]) > 0 or float(values[40]) > 0 or float(values[41]) > 0 or float(values[42]) > 0 or float(values[43]) > 0 or float(values[44]) > 0 or float(values[45]) > 0 or float(values[46]) > 0 or float(values[47]) > 0 or float(values[48]) > 0 or float(values[49]) > 0 or float(values[50]) > 0 or float(values[51]) > 0 or float(values[52]) > 0 or float(values[53]) > 0 or float(values[54]) > 0 or float(values[55]) > 0 or float(values[56]) > 0 or float(values[57]) > 0 or float(values[58]) > 0 or float(values[59]) > 0 or float(values[60]) > 0 or float(values[61]) > 0 or float(values[62]) > 0 or float(values[63]) > 0 or float(values[64]) > 0 or float(values[65]) > 0 or float(values[66]) > 0 or float(values[67]) > 0:
            outfile_not_unique_mutations.write('\t'.join([chrom, pos, ref, mut, score]) + "\t" + "\t".join(values) + "\t" +"\n")
        else:
            if float(values[1]) >= 0.08:
                outfile_unique.write('\t'.join([chrom, pos, ref, mut, score]) + "\t" + "\t".join(values) + "\t" +"\n")

            else:
                outfile_low_freq_unique.write('\t'.join([chrom, pos, ref, mut, score]) + "\t" + "\t".join(values) + "\t" +"\n")
            
            sys.stdout.write(" \r%d / %d " % (linecounter,j))
            sys.stdout.flush()


sys.stdout.write("\n")

outfile_unique.close()
outfile_low_freq_unique.close()
outfile_not_unique_mutations.close()

