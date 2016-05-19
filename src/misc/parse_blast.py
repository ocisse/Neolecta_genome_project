#parse BLAST output and write protein sequences of hits to file
sps = [] #keep species names here
with open('species_list.txt') as sps_in:
    for line in sps_in:
        sps.append(line.strip('\n'))


inf = "C:\\Users\\Anh\\Scer_IC-hemi_out" #BLAST output

import itertools as it
import csv

for f_i in range (1):

    al = [] #list containing blast output of each gene
    result = [] 
    with open(inf) as f:
        for key,group in it.groupby(f,lambda line: line.startswith('# BLASTP 2.2.28+')): #remove this key from the rest of the block
            if not key:
                group = list(group)
                al.append(group)
    for gene in al:
        print(gene[0])
        sps_pre = []#keep only 1 hit to each sps
        for i in range(4,len(gene)):
            taxid = int(gene[i][1:3])
            splt = gene[i].split('\t')
            
            if taxid not in sps_pre:
                sps_pre.append(taxid)
                result.append('>'+sps[taxid])
                with open("C:\\Users\\Anh\\proteomes\\"+str(taxid)+".csv") as proteome:
                    proteins = csv.reader(proteome,delimiter = ',')
                    for prot in proteins:
                        if splt[0] in prot[0]:
                            result.append(prot[1])
    with open(str(f_i)+'blasthits.txt','w') as outf:
        for prot in result:
            outf.write(prot)
            outf.write('\n')
