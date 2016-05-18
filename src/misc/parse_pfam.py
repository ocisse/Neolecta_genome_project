##get domains present in each prot from pfam output
path = "C:\\Users\\Anh\\neolecta\\constrained_neo\\pfamout\\"

out = []
domains = []
for i in range(1025):
    try:
        dat = []
        with open(path+str(i+1)+".pfamout") as inf:
            for line in inf:
                if "Domain annotation for each model" in line:
                    break
                dat.append(line.split())
        to_out = [dat[12][1]]
        for j in dat[13:]:
            try:
                e = float(j[0])
                if e<=0.001:
                    to_out.append(j[8])
                    to_out.append([e,' '.join(j[9:])])
            except:
                pass
    except:
        print i+1
    out.append(to_out)
    for d in range (1,len(to_out),2):
        domains.append(to_out[d])
