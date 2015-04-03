find . -type f -empty -print0 -exec rm -v "{}" \;
qsub -I -q highmem -l nodes=m01
cat Spom_FPKM_logPhase.txt | tr "\r" "\n" > x
grep -v '^>' tmp.fa | grep -o [ACTGNactg] | wc -l
qsub -d `pwd` test.sh
genewise -splice_gtag -quiet -gff -pretty -alb -hmmer /opt/cegma/2.5/data/hmm_profiles/KOG1980.hmm /tmp/genomic9523.fa
perl /rhome/ocisse/bigdata/utils/qsub-pbs.pl --maxjob 50 --interval 120 --resource walltime=100:00:00 --convert no AllJobs.sh
cat failed.txt | cut -f 2 -d ':' | cut -f 2 -d ' ' | while read line; do head $line; done
# 
~/bigdata/Project3_cema/Version1/src/fgmp/1.0/src # current loca
pushd ../~/workship/shell-novice/filesystem/nelle # go to that direct but keep track of the old 
I am in ../~/workship/shell-novice/filesystem/nelle
popd # retrun ~/bigdata/Project3_cema/Version1/src/fgmp/1.0/src

# kill force
kill -9 PID 

# useful commnands
qstatMonitor
qstat -q
qalter
qdel 34 -t 48-50

# 
./ -V

cd - # return in my old dir
pushd, popd
qsub -d direct test.sh
qsub -d `pwd` test.sh
qsub -d `pwd` -I 
echo $PBS_NODEFILE

#PBS -l nodes=1:ppn=1 

N=$PBS_ARRAYID
if [ ! $N]; then
 N=$1
fi

FILELIST=list
LINE=`head -n $FILELIST | tail -n 1`

echo "LINE $N is $LINE"


env | grep PBS


-t 1-2 run.blast-array.sh
 

/usr/bin/bp_dbsplit

#
https://github.com/sujaikumar/assemblage/blob/master/README-meloidogyne.md
