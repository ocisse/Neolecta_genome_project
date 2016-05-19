##########################################################################################################
## csv to fasta
import csv
path = "C:\\Users\\Anh\\pex\\"


infile = open(path+"seqs.csv")

data = []
reader = csv.reader(infile, delimiter = ',')
for row in reader:
    data.append(row)
infile.close()


for i in range(1,len(data)+1):
    outfile = open(path+str(i)+".txt", 'w')
    outfile.write(">"+data[i-1][0])
    outfile.write('\n')
    outfile.write(data[i-1][1])
    outfile.write('\n')
    outfile.close()


##########################################################################################################
## fasta to csv
path = "C:\\Users\\Anh\\"
infile = open(path+'48.fa')
data = []
l = None 
for line in infile:
    if line[0] == '>':
        if l:
            data.append(l)
        l = line + ' '
    else:
        l += line[:-1]
if l:
    data.append(last)

new = []
for i in data:
    if i[0]=='>':
        a = i.split('\n')
        new.append(a)


import csv
with open(path+"48.csv", "w",newline='') as fp:
    a = csv.writer(fp,delimiter = ',')
    a.writerows(new)

