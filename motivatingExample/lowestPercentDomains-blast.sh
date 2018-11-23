for i in `cat lowestPercentDomains-queries.txt`; do
    echo $i

    fastaFilename="$i.fa"
    if [ ! -s "$fastaFilename" ]; then
	grep -A 1 ">$i" singleDomains.fa > $i.fa
    fi
    
    blastOutFilename="$i.blastOut"
    if [ ! -s "$blastOutFilename" ]; then 
	cmd="blastp -db singleDomainDB -query $i.fa -evalue 1000 -num_descriptions 10000 -num_alignments 10000 > $blastOutFilename"
	echo $cmd
	eval $cmd
    fi
    
    family=`grep "^$i" relevancy.tab | cut -f 2`
    echo "family: $family"
    cmd="../benchmarkingScripts/classifyRelevance.pl  --taxon=$i  --family=$family  --blastp=$blastOutFilename  --rel=relevancy.tab  --spougeExt --spouge=$i.spouge1000  --norandomsAsIrrelevants"
    echo $cmd
    eval $cmd
done
