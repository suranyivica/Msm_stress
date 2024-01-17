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
bams_input = args[3:]
print(bams_input)
APA_index = bams_input.index("BAM_copies/APAW.rg.bam")
bams = [args[3 + int(APA_index)]]
bams_input.pop(int(APA_index))
#print(bams)

vcf_name = os.path.basename(args[1]).split("_")[0]
vcf_type_name = os.path.basename(args[1]).split("_")[0][:-1]

for bam in bams_input:
    if vcf_name in bam:
        bams.append(bam)

for bam in bams_input:
    if vcf_type_name in bam:
        if bam not in bams:
            bams.append(bam)

for bam in bams_input:
    if vcf_type_name not in bam:
        bams.append(bam)
#print(bams)


if len(args) == 0:
    sys.exit("Usage: python3 getaf.py [position file] [output file for all mutations] [bam files]")
if not all(x.endswith("bam") for x in bams):
    sys.exit("Usage: python3 getaf.py [position list] [output file for all mutations] [bam files]" + " ".join([bams]))

j = int(subprocess.check_output(['wc', '-l', args[1]]).split()[0])

outfile = open(args[2], "w")
bam_h = list(map(os.path.basename, bams))
bam_names = []
for elements in bam_h:
    names_parts = elements.split(".")
    bam_names.append(names_parts[0])
outfile.write("\t".join(["chromosome", "pos", "ref", "mut"]) + "\t" + "\t".join(bam_names) + "\n")

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

        outfile.write('\t'.join([chrom, pos, ref, mut]) + "\t" + "\t".join(values) + "\t" +"\n")
            #sys.stdout.write(" \r%d / %d " % (linecounter,j))
            #sys.stdout.flush()



sys.stdout.write("\n")

outfile.close()
