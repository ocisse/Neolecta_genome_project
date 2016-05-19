import csv
from subprocess import call
from Bio import AlignIO
from Bio import Phylo

pezizo_taxids = [11,9,21,5,13,35,18,26,1,14,2,4]

sps = [] #keep species names here
with open('species_list.txt') as sps_in:
    for line in sps_in:
        sps.append(line.strip('\n'))

path1 = "C:\\Users\\Anh\\dyn+pex\\pex\\hmm\\"


##get seq IDs
all_ids = [] #stores the ID of best hit
all_es = [] # stores e-val of best hit
out = []

for i in range(1,430):
    ids = ['-' for t in range(43)]
    es = ['-' for g in range(43)]
    inf = open(path1+str(i)+'.ancienthmmout')
    read = []
    for line in inf:
        if '------ inclusion threshold ------' in line:
            break
        else:
            read.append([x for x in line.split(' ') if x!=''])     
    inf.close()
    sps_pre = []#keep only 1 hit to each sps
    for entry in read[14:]:
        if len(sps_pre)>40:
            break
        
        elif float(entry[0])<=0.0001:
            try:
                taxid = int(entry[8][1:3])

                if taxid not in sps_pre:
                    ids[taxid] = entry[8]
                    es[taxid] = entry[0]
                    sps_pre.append(taxid)
            except:
                continue
                
    all_ids.append(ids)
    all_es.append(es)
          
####remove homologs that fails the reciprocal best hit test
for column in range(43):
    col = []
    for row in range(len(all_ids)):
        cur = all_ids[row][column]
        if cur == '-':
            col.append(cur)
        else:
            if cur not in col:
                col.append(cur)
            else:
                remove_this = col.index(cur)
                if float(all_es[row][column])<float(all_es[remove_this][column]):
                    all_ids[remove_this][column]='-'
                    col[remove_this]='-'
                    col.append(cur)
                elif float(all_es[row][column])>float(all_es[remove_this][column]):
                    all_ids[row][column]='-'
                    col.append('-')

##write                
for r in range(len(all_ids)):
    allseqs = []
    for c in range(len(all_ids[0])):
        if all_ids[r][c]!='-':
            taxid = int(all_ids[r][c][1:3])
            with open("C:\\Users\\Anh\\proteomes\\"+str(taxid)+".csv") as proteome:
                proteins = csv.reader(proteome,delimiter = ',')
                collect_seqs = []
                for prot in proteins:
                    if all_ids[r][c] in prot[0]:
                        collect_seqs.append(prot[1])   ##have to do this because bloody Nirr is not nr
                if len(collect_seqs)==1:
                    allseqs.append('>'+sps[taxid])
                    allseqs.append(collect_seqs[0])

                elif len(collect_seqs)>1:
                    s = collect_seqs[0]
                    for seqs in collect_seqs[1:]:
                        if len (seqs)>len(s):
                            s = seqs
                    allseqs.append('>'+sps[taxid])
                    allseqs.append(s)


    with open(path1+str(r+1)+'.seqs','w') as outf:
        for d in allseqs:
            outf.write(d)
            outf.write('\n')
            
            
##get MSA, trim, build the tree with PhyML then come back here for the following steps


##extract branch lengths/subs rates
dist = [sps+['pez_avg']]

pez_10char = [] #10-character names of pezizo species - coz .phylip files accept only <= 10 char for seq identifiers
for ind in pezizo_taxids:
    pez_10char.append(sps[ind][:10])


for i in range(1,14):
    tree = Phylo.read(path1+str(i)+".phy_phyml_tree",'newick')
    
    leaves = tree.get_terminals()
    pez_leaves = []  #pez species present in this particular tree. need this because not all trees contain the same pez species
    for leave in leaves:
        for sp in pez_10char:
            if leave.name==sp:
                pez_leaves.append(sp)
                break
    pezizo = tree.common_ancestor({"name":spsname} for spsname in pez_leaves)
    
    to_dist = ['-' for i in range(43)]
    for leave in leaves:
        for s in sps:
            if leave.name==s:#[:10]:
                to_dist[sps.index(s)]=tree.distance(pezizo,leave.name)#/tree.distance(pezizo, "N.irregula")
                break

    pez = 0 #accumulates total branch lengths of all pezizo sps to the pez root, to calculate the average later, as a proxy of the 'radius' of this cluster
    for sp in pez_leaves:
        pez+=tree.distance(sp, pezizo)

    
    pez_avg = pez/len(pez_leaves)
    to_dist.append(pez_avg)
    dist.append(to_dist)
##
##
order = [42,22, 11, 9,4,5,21,13,2,35,14,1,18,26,24,6,32,34,33,28,12,20,15,25,30, #rearrange the output from aphabetical order of sps names to phylogenetic rel-based order
         16,29,31,27,23,17,10,8,7,19,0,3,
         36, 37, 38, 39, 40,41]

rearrange = []
for ent in dist:
    to_rearrange = []
    for number in order:
        to_rearrange.append(ent[number])
    to_rearrange.extend(ent[-1:])
    rearrange.append(to_rearrange)

only_necessary_cols_ids = [0,9,13,14,15,20,21,23,25,30,27,29,31,33,34,36,38,42,40,37,41,44]

to_csv = []
for ent in rearrange:
    t = []
    for number in only_necessary_cols_ids:
        t.append(ent[number])
    to_csv.append(t)


outf = open(path1+'distmuscle.csv','w',newline='')
write = csv.writer(outf, delimiter = ',')
write.writerows(to_csv)
outf.close()
