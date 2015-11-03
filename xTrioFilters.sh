#!/bin/bash
#$ -cwd -l mem=1G,time=3:: -N FilterTrio

usage="xTrioFilters.sh -v <vcf file> -t <Trio Table> -n <Output Directory> -a <The line of the Trio Table to analyse> -p <Additional Parameters> -D <ignore de novo quality> -H <this message>"

BadDeN=false
ArrNum=1
#get arguments
while getopts v:t:n:a:p:DH opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) TrioFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        a) ArrNum="$OPTARG";;
        p) AddPrm="$OPTARG";;
        D) BadDeN="true";;
        H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`
TrioFil=`readlink -f $TrioFil`


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
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.denovo --het $Proband --ref $Father,$Mother"
if [[ "$BadDeN" == "false" ]]; then CMD=$CMD" -D"; fi #otherwise de novos will be run with default filters
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.denovo.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.denovo.tsv &
fi

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.AR --alt $Proband --het $Father,$Mother"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.AR.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.AR.tsv &
fi

#X linked - male proband
echo "X linked - male proband.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.X-linked --alt $Proband --het $Mother --ref $Father -X"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.X-linked.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.X-linked.tsv &
fi

#Autosomal Dominant - paternal inheritance
echo "Autosomal Dominant - paternal inheritance.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.AD-paternal  --het $Proband,$Father --ref $Mother"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.AD-paternal.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.AD-paternal.tsv &
fi

#Autosomal Dominant - maternal inheritance
echo "Autosomal Dominant - maternal inheritance.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.AD-maternal  --het $Proband,$Mother --ref $Father"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.AD-maternal.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.AD-maternal.tsv &
fi

#compound heterozygous
echo "Compund heterozygous.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.tempheppat  --het $Proband,$Father --ref $Mother -P  -f 0.03"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.Trio.temphepmat  --het $Proband,$Mother --ref $Father -P  -f 0.03"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$FamNam.Trio.temphepmat.tsv")
pathet <- read.delim("$FamNam.Trio.tempheppat.tsv")

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$FamNam.Trio.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
cat $FamNam.Trio.tempheppat.log $FamNam.Trio.temphepmat.log > $FamNam.Trio.compound_heterozygous.log
LEN=`cat $FamNam.Trio.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.compound_heterozygous.tsv &
fi


rm -rf $FamNam.Trio.temp*
