## input: list of genes with their presence/absence in each species
## output: the same list of genes clustered by their presence/absence

import csv
import numpy as np
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt
path = "C:\\Users\\Anh\\transport_routes_2nd_hmm\\sto\\phy\\"
cluster_in = path + "cluster_in.csv"
cluster_out = path + "cluster_out.csv"
PA_matrix = path + "dist_trimmed.csv"


##clustering
inf = open(cluster_in)
read = csv.reader(inf, delimiter = ',')
col = []
for line in read:
    col.append(line)
inf.close()
ar = [k[:-2] for k in col[1:]]
arr = np.array(ar)
lkg = hac.linkage(arr,'complete')
hac.dendrogram(lkg)
##plt.show()

##copy the presence/absence indicators (1/0) of each gene to the new position in the cluster
order = []
for i in hac.leaves_list(lkg):
    order.append(i)

readPA = open(PA_matrix)
PAreader = []
r = csv.reader(readPA,delimiter=',')
for line in r:
    PAreader.append(line)
readPA.close()

##write
out = []
out.append(PAreader[0])
for index in order:
    out.append(PAreader[index+1])
	
with open(cluster_out,'w',newline='') as outfile:
    writer = csv.writer(outfile, delimiter = ',')
    writer.writerows(out)
