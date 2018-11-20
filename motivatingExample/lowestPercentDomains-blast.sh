for i in `cat lowestPercentDomains.txt`; do
    echo $i
    grep -A 1 singleDomain.fa > $i.fa
    cmd="blastp -db singleDomainDB -query $i.fa -evalue 1000 -num_descriptions 10000 -num_alignments 10000 > $i.blastOut"
    echo $cmd
    eval $cmd
    family=`grep "^$i" | cut -f 2`
    echo "family: $family"
    cmd="../benchmarkingScripts/classifyRelevance.pl  --taxon=$i  --family=$family  --blastp=$i.blastOut  --rel=relevancy.tab  --spougeExt=$i.spouge1000"
    echo $cmd
    eval $cmd
done