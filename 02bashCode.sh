#! /usr/bin/bash
#
#	github.com/pjbiggs
#
# 	last edited 2023-09-14 
# 
#	Code to run a short conversion for the sequences from dada2 into a text file for input back in the R environment.
#	This is run from the same folder where the dada2 results are.
# 	It requires that cd-hit be accessible in the path, and so to are all the extra scripts for working with cd-hit output

 
grep -v '^x' dada_seqs.txt | awk '{print ">seq"$1"\n"$2}' > dada_seqs.fa
grep -v '^x' dada_seqs.txt | awk '{print "seq"$1"\t"$2}' > dada_seqsMod.txt
cd-hit -i dada_seqs.fa -o dada_clust.fa
head -n 10 dada_seqs.fa
perl clstr_sql_tbl.pl dada_clust.fa.clstr dada_clust.txt
cat dada_clust.txt | perl -lpe 's/Hom/C\. \hominis/g;s/Par/C\. \parvum/g' > species.txt

