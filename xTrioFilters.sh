#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterTrio

usage="(-t 1-X) xTrioFilters.sh -v <vcf file> -t <Trio Table> -n <Output Directory> -p <Additional Parameters> -D <ignore de novo quality> -H <this message>"

BadDeN=false
#get arguments
while getopts v:t:n:p:DH opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) TrioFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        p) AddPrm="$OPTARG";;
        D) BadDeN="true";;
        H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
TrioFil=`readlink -f $TrioFil`


if [[ "$SGE_TASK_ID" != "undefined" ]]; then
        ArrNum=$SGE_TASK_ID
else
        ArrNum=1
fi

echo $ArrNum

FamNam=`cut -f 1 $TrioFil | head -n $ArrNum | tail -n 1`
Proband=`cut -f 2 $TrioFil | head -n $ArrNum | tail -n 1`
Father=`cut -f 3 $TrioFil | head -n $ArrNum | tail -n 1`
Mother=`cut -f 4 $TrioFil | head -n $ArrNum | tail -n 1`
Extras=`cut -f 5 $TrioFil | head -n $ArrNum | tail -n 1`


echo $FamNam
if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#De novo
echo "De Novo Filtering.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.denovo --het $Proband --ref $Father,$Mother"
if [[ "$BadDeN" == "false" ]]; then CMD=$CMD" -D"; fi #otherwise de novos will be run with default filters
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.denovo.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.denovo.tsv
fi

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AR --alt $Proband --het $Father,$Mother"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.AR.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.AR.tsv
fi

#X linked - male proband
echo "X linked - male proband.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.X-linked --alt $Proband --het $Mother --ref $Father -X"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.X-linked.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.X-linked.tsv
fi

#Autosomal Dominant - paternal inheritance
echo "Autosomal Dominant - paternal inheritance.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AD-paternal  --het $Proband,$Father --ref $Mother"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.AD-paternal.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.AD-paternal.tsv
fi

#Autosomal Dominant - maternal inheritance
echo "Autosomal Dominant - maternal inheritance.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AD-maternal  --het $Proband,$Mother --ref $Father"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.AD-maternal.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.AD-maternal.tsv
fi

#compound heterozygous
echo "Compund heterozygous.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.tempheppat  --het $Proband,$Father --ref $Mother"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.temphepmat  --het $Proband,$Mother --ref $Father"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$FamNam.temphepmat.tsv")
pathet <- read.delim("$FamNam.tempheppat.tsv")

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$FamNam.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
cat $FamNam.tempheppat.log $FamNam.temphepmat.log > $FamNam.filter.compound_heterozygous.log
LEN=`cat $FamNam.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.compound_heterozygous.tsv
fi


rm -rf *temp*
