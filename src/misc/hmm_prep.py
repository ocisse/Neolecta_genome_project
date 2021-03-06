inseqpath = "C:\\Users\\Anh\\neolecta\\dyn+pex\\"



####get MSA
from subprocess import call
for i in range(1050):
    muscle = "C:\\Users\\Anh\\PhyloKits\\muscle.exe -in "
    inp = inseqpath+ str(i)+"seqs.txt -out "
    outp= inseqpath + str(i)+ "aln"

    command = muscle+inp+outp
    call(command)


######trim MSA
##from subprocess import call
for i in range(1050):
    trimal = "C:\\Users\\Anh\\PhyloKits\\trimAl\\bin\\trimal -in "
    inp = inseqpath + str(i)+ 'aln -out '
    outp = inseqpath+ str(i) + 'alntrim.txt -gt 0.8 -cons 0.5'
    command = trimal + inp + outp
    call(command)

#######convert fasta to stockholm
from Bio import AlignIO
for i in range(1050):
    inf = open(inseqpath+ str(i)+'alntrim.txt')
    alignment = AlignIO.parse(inf, "fasta")
    outf = open(inseqpath+ str(i)+'.sto', 'w')
    AlignIO.write(alignment, outf, "stockholm")
    inf.close()
    outf.close()
