from __future__ import division
from subprocess import call
from scipy.spatial import distance
from scipy import signal
import numpy as np
import pylab as plt
import math

pez = ['>N.crassa','>A.oligospora','>B.graminis','>F.oxysporum','>M.oryzae','>N.fischeri','>P.confluens','>S.sclerotiorum','>A.nidulans']

path = "C:\\Users\\Anh\\plots\\PEX\\"


#read substitution matrix
read = []
inf = open('C:\\Users\\Anh\\BLOSUM62.txt')
for line in inf:
    read.append(line.strip('\n').split('\t'))
inf.close()
row  = read[0]
col  = [s[0] for s in read]


#look up the score for substitution of c2 by c1
def get_score(c1,c2):
    if c1=='-' and c2=='-':
        return 0
    elif c1=='-' or c2=='-':
        return -3 #gap penalty = 3
    return int(read[row.index(c1)][col.index(c2)])


def cal_dst(pez,other): #pez is a list of pez char at a particular position (in a column). same for other
    #this function returns the combined scores of the substitution of each Pezizo character by another character in the 'other' group
    score = 0
    for c_p in pez:
        for c_o in other:
            score+=get_score(c_p,c_o)
    return score/len(other)

def normalize(lst,maxi,mini):
    out = []
    for num in lst:
        abs_max = maxi if maxi>abs(mini) else abs(mini)
        norm = abs(num)/abs_max
        if num>=0:
            out.append(norm)
        else:
            out.append(0-norm)
    return out

#check whether the column is conserved
def consd(t1,t2,t3):
    l = len(t1)+len(t2)+len(t3)
    c = [0 for i in range(7)]
    all_c = t1
    all_c.extend(t2)
    all_c.extend(t3)
    for char in all_c:
        for group in range(len(similar)):
            if char in similar[group]:
                c[group] +=1
                break
    return True if max(c)>0.7*l else False

#read input alignment
def read_input(filename):
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

def split_lst(xlst,ylst):
    s = []
    x,y = [],[]
    for i in range(len(xlst)-1):
        if xlst[i]!=xlst[i+1]-1:
            s.append(i)

    x.append(xlst[0:s[0]+1])
    for i in range(len(s)-1):
        x.append(xlst[s[i]+1:s[i+1]+1])
    x.append(xlst[s[-1]+1:])

    for l in x:
        c = []
        for b in l:
            c.append(ylst[b])
        y.append(c)
    return x,y


peroxins = [] ##insert names of seqs here

raw_scores = []
dat = [] # stores lists of y_ scores after SG filter, before normalization
gap_dat = []
skip = []
conserved = []
conserved_in_CM = []

for idx in range(1,19):
    print(idx)
    collect = read_input(path+str(idx)+'aln')
    cons = []
    CM_cons = []
    scores = [] #store the raw scores obtained using cal_dst here
    gap = []
    gap_both = []
    ncra_seq = ''
    pez_seqs,yeast_seqs,neo_seq = [],[],[]
    for s in collect:
        if s[0]in pez:
            pez_seqs.append(s[1])
        elif s[0]=='>N.irregularis':
            neo_seq.append(s[1])
        else:
            yeast_seqs.append(s[1])

    current_actual_pos = 0 #with respect to N.irr seq
    for pos in  range(len(ncra_seq)):
        if neo_seq[0][pos]!='-':

            pc = [s[pos] for s in pez_seqs]
            yc = [s[pos] for s in yeast_seqs]
            nc = [s[pos] for s in neo_seq]

            p_gap = True if pc.count('-')>=math.ceil(0.5*len(pc)) else False
            h_gap = True if yc.count('-')>=math.ceil(0.5*len(yc)) else False
            n_gap = True if nc[0]=='-' else False

            ys = cal_dst(nc,yc)
            ns = cal_dst(nc,pc)
            sc = ns-ys
            scores.append(sc)
##            print yc,pc,sc
            if consd(pc,nc,yc):
                cons.append(current_actual_pos)
            if CMconsd(nc,pc):
                CM_cons.append(current_actual_pos)
                
            if n_gap:
                gap_both.append(current_actual_pos)
            if not n_gap and (h_gap or n_gap):
                gap.append(current_actual_pos)
                
            current_actual_pos += 1
    dat.append(signal.savgol_filter(scores,41,5,mode='mirror'))
    gap_dat.append(gap)
    skip.append(gap_both)
    raw_scores.append(scores)
    conserved.append(cons)
    conserved_in_CM.append(CM_cons)
    

maxis = [max(l) for l in raw_scores]
minis = [min(m) for m in raw_scores]
global_max = max(maxis)
global_min = min(minis)
correct_pos = 1.67/global_max
correct_neg = 1.67/global_max

for indx in range(len(raw_scores)):
    print(peroxins[indx])
    
    raw_norm = normalize(raw_scores[indx],global_max,global_min)
    avg_norm = normalize(dat[indx],global_max,global_min)
    
    for i in range(len(raw_scores[indx])):
        corrected = max(0,raw_norm[i]-correct_pos) if raw_norm[i]>0 else min(0,raw_norm[i]+correct_neg)
        if i in conserved[indx]:
            print(i,raw_norm[i],avg_norm[i],'. ',1)

        else:
            if i in skip[indx]:
                print(i,'.','.','.','.')
            else:
                if i in gap_dat[indx]:
                    print(i,0,avg_norm[i],raw_norm[i],'.')

                else:
                    print(i,raw_norm[i],avg_norm[i],'.','.')
