import csv
sps = [] #keep species names here
with open('species_list.txt') as sps_in:
    for line in sps_in:
        sps.append(line.strip('\n'))

         
path = "C:\\Users\\Anh\\hmmoutput\\"

##get seqs

for i in range(8,12):
    
    tofile = []
    inf = open(path+str(i)+'.hmmout')
    read = []
    for line in inf:
        if '------ inclusion threshold ------' in line:
            break
        else:
            read.append([x for x in line.split(' ') if x!=''])     
    inf.close()
    sps_pre = []#keep only 1 hit to each sps
    for entry in read[14:]:
        if float(entry[0])<=0.001:
            taxid = int(entry[8][1:3])
            if taxid not in sps_pre:
                with open("C:\\Users\\Anh\\proteomes\\"+str(taxid)+".csv") as proteome:
                    proteins = csv.reader(proteome,delimiter = ',')
                    for prot in proteins:
                        if entry[8] in prot[0]:
                            tofile.append(['>'+sps[taxid], prot[1]])                                          
                sps_pre.append(taxid)

    with open(path+str(i)+'seqs.txt','w') as outf:
        for d in tofile:
            outf.write(d[0])
            outf.write('\n')
            outf.write(d[1])
            outf.write('\n')
            
            

####get MSA from muscle
inseqpath = "C:\\Users\\Anh\\seqs\\"

from subprocess import call
from Bio import AlignIO
from Bio import Phylo
for i in range(1050):
    muscle = "C:\\Users\\Anh\\PhyloKits\\muscle -in "
    inp = inseqpath+ "hmmseqs"+ str(i)+'seqs.txt -out '
    outp= inseqpath + str(i)+"aln"

    command = muscle+inp+outp
    call(command)

#convert fasta to phylip
for i in range(1050):
    inf = open(inseqpath+str(i)+'aln')
    alignment = AlignIO.parse(inf, "fasta")
    outf = open(inseqpath+str(i)+'.phy', 'w')
    AlignIO.write(alignment, outf, "phylip")
    inf.close()
    outf.close()
