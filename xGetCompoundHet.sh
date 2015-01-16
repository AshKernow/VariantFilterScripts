#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N GetCompHet

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts p:m:o: opt; do
    case "$opt" in
        p) PatHet="$OPTARG";;
        m) MatHet="$OPTARG";;
        o) OutNam="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$PatHet")
pathet <- read.delim("$MatHet")

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],match(colnames(mathet), colnames(pathet))]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$OutNam.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT

echo "Paternal Filtering" > $OutNam.compound_heterozygous.log
echo "---------------------------------------------------" >> $OutNam.compound_heterozygous.log
cat ${PatHet/tsv/log} >> $OutNam.compound_heterozygous.log
echo "" >> $OutNam.compound_heterozygous.log
echo "===================================================" >> $OutNam.compound_heterozygous.log
echo "" >> $OutNam.compound_heterozygous.log

echo "Maternal Filtering" >> $OutNam.compound_heterozygous.log
echo "---------------------------------------------------" >> $OutNam.compound_heterozygous.log
cat ${MatHet/tsv/log} >> $OutNam.compound_heterozygous.log
echo "" >> $OutNam.compound_heterozygous.log
echo "===================================================" >> $OutNam.compound_heterozygous.log
echo "" >> $OutNam.compound_heterozygous.log
OvLen=`cut -f 6  $OutNam.compound_heterozygous.tsv |  sort | uniq | wc -l`
OvLen=$(( OvLen - 1 ))
echo "Number of genes with possible compound heterozyogus variants: $OvLen" >> $OutNam.compound_heterozygous.log
OvLen=`cat $OutNam.compound_heterozygous.tsv | wc -l`
OvLen=$(( OvLen - 1 ))
echo "Number of variants: $OvLen" >> $OutNam.compound_heterozygous.log
