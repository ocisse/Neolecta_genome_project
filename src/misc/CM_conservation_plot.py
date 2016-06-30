##This script was used to generate the plots in Figs. 5c and 6c.
##input : a multiple sequence alignment (MSA) consisting of homologs from clade A, clade B, and sequence x from a different clade, acting as a reference.
##        In the Neolecta paper, A = Pezizomycotina, B = Saccharomycotina, x is from N. irregularis.
##output: a similarity score for each amino acid position in x. A positive score means that position in x is more similar to 
##        corresponding amino acids in group A in the MSA. A negative score means that position in x is more similar to corresponding
##        amino acids in group B in the MSA.
##
##        To smoothen the curve/reduce noise, the Savitzky-Golay filter is applied with window size = 41, polynomial order 5.
##        The scores will then be plotted against amino acid positions in x, with the dotted line representing regions in x that are gaps
##        in >=70% of the sequences in the input MSA.


from __future__ import division
from subprocess import call
from scipy.spatial import distance
from scipy import signal
import numpy as np
import pylab as plt
import math

##names of species in clade A (Pezizomycotina) in the MSA
pez = ['>N.crassa','>A.oligospora','>B.graminis','>F.oxysporum','>M.oryzae','>N.fischeri','>P.confluens','>S.sclerotiorum','>A.nidulans']

##names of the protein sequences of interest. used to name output files
proteins = []

##directory of input MSAs. Output will be written here.
path = "C:\\Users\\Anh\\plots\\PEX\\"

##similar amino acids, grouped by their biochemical properties. Used to determine overall conservation in each column of the MSA
##used to determine conserved positions (positions with a red line in the background in Fig 5,6)
similar = [['G','A','V','L','I'],
 ['F','Y','W'],
 ['C','M'],
 ['S','T'],
 ['K','R','H'],
 ['D','E','N','Q'],
 ['P']]

##read the substitution matrix
mat = []
inf = open('C:\\Users\\Anh\\plots\\BLOSUM62.txt')
for line in inf:
    mat.append(line.strip('\n').split('\t'))
inf.close()
row  = mat[0]
col  = [s[0] for s in mat]

##########################################################################################################################################################
## Helper functions                                                                                                                                     ##
##########################################################################################################################################################

def get_score(c1,c2):
##    input: 2 characters c1 and c2
##    output: score for substituting c2 by c1
    if c1=='-' and c2=='-': ## if both characters are gaps, return 0
        return 0
    elif c1=='-' or c2=='-': ## apply gap penalty = 3 if only 1 character is a gap
        return -3
    return int(read[row.index(c1)][col.index(c2)]) ##otherwise, look for the score of subtituting c2 by c1 in the matrix 'mat'


def cal_dst(A,B):
##    input: 2 lists, A and B. A is a list of char at a particular position (in a column) in clade A. Likewise for B.
##    output: average sum-of-pairs score divided by the number of sequences in each clade A and B (this provides normalization in case
##    the numbers of sequences are different in the 2 input clades)  
    score = 0
    for c_p in pez:
        for c_o in other:
            score+=get_score(c_p,c_o)
    return score/len(other)


def normalize(lst,maxi,mini):
##    input: a list of numbers (lst), a global max and a global min (these must not be inside the range of lst)
##    output: normalized version of the original list, based on global max and min. the max absolute value of global max and min will be used
##    as the max and -min of the returned list
    out = []
    for num in lst:
        abs_max = maxi if maxi>abs(mini) else abs(mini)
        norm = abs(num)/abs_max
        if num>=0:
            out.append(norm)
        else:
            out.append(0-norm)
    return out


def consd(t1,t2,t3):
##    input: 3 lists of characters in the same alignment column. 3 lists because the species are originally divided in 3 groups.
##    output: True if at more than 70% of the characters belong to the same biochemical group. False otherwise
   
    l = len(t1)+len(t2)+len(t3)
    c = [0 for i in range(7)] # 7 biochemical groups (refer to 'similar' above)
    all_c = t1
    all_c.extend(t2)
    all_c.extend(t3)
    for char in all_c:
        for group in range(len(similar)):
            if char in similar[group]:
                c[group] +=1
                break
    return True if max(c)>0.7*l else False


def read_input(filename):
##    input: an MSA in fasta format
##    output: a list of pairs [species_name, sequence as it appears in the MSA]
    
    infile = open(filename)
    data = []
    last = None
    for line in infile:
        if line[0] == '>':
            if last:
                data.append(last)
            last = line + ' '
        else:
            last += line[:-1]
    if last:
        data.append(last)
    coll = []
    for li in data:
        a = li.split('\n ')
        coll.append(a)
    return coll

##########################################################################################################################################################
## main function                                                                                                                                        ##
##########################################################################################################################################################

raw_scores = [] # list of lists. A child list stores the raw score returned by cal_dst(A,B) for each alignment position in the MSA of a protein
                #ignore positions that are gap in the reference sequence (x)
dat = [] # stores lists of y_ scores after Savitzky-Golay filter, before normalization
gap_dat = []
skip = []
conserved = []

for idx in range(len(proteins)):
    collect = read_input(path+str(idx)+'aln') #read input alignment
    cons = []
    scores = [] #store the raw scores obtained using cal_dst here
    gap = []
    gap_both = []
    ncra_seq = ''
    pez_seqs,yeast_seqs,neo_seq = [],[],[]
    for s in collect:
    #assign the sequences to their groups: Pezizomycotina, Saccharomycotina, and Neolecta
        if s[0]in pez:
            pez_seqs.append(s[1])
        elif s[0]=='>N.irregularis':
            neo_seq.append(s[1])
        else:
            yeast_seqs.append(s[1])

    current_actual_pos = 0 #with respect to Neoleta seq (ignore gaps, hence 'actual')
    for pos in  range(len(neo_seq)):
        if neo_seq[0][pos]!='-': #only consider positions that are not gaps in Neolecta
            #put the characters in that alignment column into lists corresponding to the 3 groups
            pc = [s[pos] for s in pez_seqs]
            yc = [s[pos] for s in yeast_seqs]
            nc = [s[pos] for s in neo_seq]

            #check whether the position is a gap  in either group (if >50% are gaps in the seqs of a group, the postition is considered a gap in that group)
            p_gap = True if pc.count('-')>=math.ceil(0.5*len(pc)) else False
            h_gap = True if yc.count('-')>=math.ceil(0.5*len(yc)) else False

            #calculate the difference between Neolecta-Pezizomycotina and Neolecta-Saccharomycotina sums-of pairs. Append that to 'scores'
            ys = cal_dst(nc,yc)
            ns = cal_dst(nc,pc)
            sc = ns-ys
            scores.append(sc)

            #check if that position is conserved in all 3 groups. if yes, add the postion number to 'cons' (for the red background)
            if consd(pc,nc,yc):
                cons.append(current_actual_pos)

            #if the position is not a gap in Neolecta, and a gap in either of the other 2 clades, add the position number to 'gap' (for the dotted line)
            if not n_gap and (h_gap or n_gap):
                gap.append(current_actual_pos)
            
            #move on to the next position in the alignment   
            current_actual_pos += 1

##  after each protein is processed, add its scores after smoothing with Savitzky-Golay filter to 'dat',
##  add the record of postitions that are gap in the other clades to 'gap_dat', add the record of conserved positions to 'conserved'
    dat.append(signal.savgol_filter(scores,41,5,mode='mirror'))
    gap_dat.append(gap)
    conserved.append(cons)
    
##  to normalize the level of conservation across all the protein MSA in the analysis, find the max among the max scores (global_max)
##  from each protein MSA. Likewise for the min.
maxis = [max(l) for l in raw_scores]
minis = [min(m) for m in raw_scores]
global_max = max(maxis)
global_min = min(minis)

##########################################################################################################################################################
## print the output                                                                                                                                     ##
##########################################################################################################################################################

##format of printed output:
##protein name
## if the position is not a gap:
##  position number         score       '.'         1 if conserved, else '.'
##if the position is a gap:
##  position number         '.'           score       '.'

for indx in range(len(raw_scores)):
    print(proteins[indx])    
    avg_norm = normalize(dat[indx],global_max,global_min)    
    for i in range(len(dat[indx])):
        if i in conserved[indx]:
            print(i,avg_norm[i],'. ',1)
        else:
            if i in gap_dat[indx]:
                print(i,'.',avg_norm[i],'.')
            else:
                print(i,avg_norm[i],'.','.')


##########################################################################################################################################################
## calculate the score of randomly generated alignment column                                                                                           ##
##########################################################################################################################################################

##randomly generate 1,000,000 alignment columns based on the amino acid frequencies specified by BLOSUM62
##calculate the score for each column
##write all the scores to an output file and print the values of 2.5th and 97.5th percentile

aa = ['M','I','L','V','F','C','A','G','P','T','S','Y','W','Q','N','H','E','D','K','R']

##BLOSUM62 background frequencies
freq =[0.025,0.068,0.099,0.073,0.047,0.025,0.074,0.074,0.039,0.051,0.057,0.034,0.013,0.034,0.045,0.026,0.054,0.054,0.058,0.052]

cumulative_freq = []
counts = [0 for c in aa]
running_total = 0
results = [] #stores the scores of randomly generated columns

##write the scores 
outf = open("C:\\Users\\Anh\\background.txt",'w')

def get_random_char():
    n = random.random()
    c = ''
    for j in range(0,len(cumulative_freq)):
        if n<cumulative_freq[j]:
            c = aa[j]
            break
    return c,j

for f in freq:
    running_total += f
    cumulative_freq.append(running_total)


for i in range(1000000):
    neo = []
    pez = []
    hem = []
    n,num = get_random_char()
    neo.append(n)
    counts[num] +=1
    for j in range(9):
        p,no = get_random_char()
        pez.append(p)
        counts[no] +=1
        h,no = get_random_char()
        hem.append(h)
        counts[no] +=1
    r = cal_dst(neo,pez)-cal_dst(neo,hem)
    results.append(r)
    outf.write(str(r))

outf.close()

arr = np.array(results)
print ('2.5th percentile: '+ str(np.percentile(arr, 2.5)))
print ('2.5th percentile: '+ str(np.percentile(arr, 97.5)))

##BLOSUM62.txt
##	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V	B	Z	X	*
##A	4	-1	-2	-2	0	-1	-1	0	-2	-1	-1	-1	-1	-2	-1	1	0	-3	-2	0	-2	-1	0	-4
##R	-1	5	0	-2	-3	1	0	-2	0	-3	-2	2	-1	-3	-2	-1	-1	-3	-2	-3	-1	0	-1	-4
##N	-2	0	6	1	-3	0	0	0	1	-3	-3	0	-2	-3	-2	1	0	-4	-2	-3	3	0	-1	-4
##D	-2	-2	1	6	-3	0	2	-1	-1	-3	-4	-1	-3	-3	-1	0	-1	-4	-3	-3	4	1	-1	-4
##C	0	-3	-3	-3	9	-3	-4	-3	-3	-1	-1	-3	-1	-2	-3	-1	-1	-2	-2	-1	-3	-3	-2	-4
##Q	-1	1	0	0	-3	5	2	-2	0	-3	-2	1	0	-3	-1	0	-1	-2	-1	-2	0	3	-1	-4
##E	-1	0	0	2	-4	2	5	-2	0	-3	-3	1	-2	-3	-1	0	-1	-3	-2	-2	1	4	-1	-4
##G	0	-2	0	-1	-3	-2	-2	6	-2	-4	-4	-2	-3	-3	-2	0	-2	-2	-3	-3	-1	-2	-1	-4
##H	-2	0	1	-1	-3	0	0	-2	8	-3	-3	-1	-2	-1	-2	-1	-2	-2	2	-3	0	0	-1	-4
##I	-1	-3	-3	-3	-1	-3	-3	-4	-3	4	2	-3	1	0	-3	-2	-1	-3	-1	3	-3	-3	-1	-4
##L	-1	-2	-3	-4	-1	-2	-3	-4	-3	2	4	-2	2	0	-3	-2	-1	-2	-1	1	-4	-3	-1	-4
##K	-1	2	0	-1	-3	1	1	-2	-1	-3	-2	5	-1	-3	-1	0	-1	-3	-2	-2	0	1	-1	-4
##M	-1	-1	-2	-3	-1	0	-2	-3	-2	1	2	-1	5	0	-2	-1	-1	-1	-1	1	-3	-1	-1	-4
##F	-2	-3	-3	-3	-2	-3	-3	-3	-1	0	0	-3	0	6	-4	-2	-2	1	3	-1	-3	-3	-1	-4
##P	-1	-2	-2	-1	-3	-1	-1	-2	-2	-3	-3	-1	-2	-4	7	-1	-1	-4	-3	-2	-2	-1	-2	-4
##S	1	-1	1	0	-1	0	0	0	-1	-2	-2	0	-1	-2	-1	4	1	-3	-2	-2	0	0	0	-4
##T	0	-1	0	-1	-1	-1	-1	-2	-2	-1	-1	-1	-1	-2	-1	1	5	-2	-2	0	-1	-1	0	-4
##W	-3	-3	-4	-4	-2	-2	-3	-2	-2	-3	-2	-3	-1	1	-4	-3	-2	11	2	-3	-4	-3	-2	-4
##Y	-2	-2	-2	-3	-2	-1	-2	-3	2	-1	-1	-2	-1	3	-3	-2	-2	2	7	-1	-3	-2	-1	-4
##V	0	-3	-3	-3	-1	-2	-2	-3	-3	3	1	-2	1	-1	-2	-2	0	-3	-1	4	-3	-2	-1	-4
##B	-2	-1	3	4	-3	0	1	-1	0	-3	-4	0	-3	-3	-2	0	-1	-4	-3	-3	4	1	-1	-4
##Z	-1	0	0	1	-3	3	4	-2	0	-3	-3	1	-1	-3	-1	0	-1	-3	-2	-2	1	4	-1	-4
##X	0	-1	-1	-1	-2	-1	-1	-1	-1	-1	-1	-1	-1	-1	-2	0	0	-2	-1	-1	-1	-1	-1	-4
##*	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	1
