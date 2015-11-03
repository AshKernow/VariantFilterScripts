#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterPC

ArrNum=1
usage="xProbandParentFilters.sh -v <vcf file> -t <Filter Table> -n <Output Directory> -a <The line of the Trio Table to analyse> -p <Additional Parameters> -H <this message>"

#get arguments
while getopts v:t:n:a:p:H opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) PartFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        a) ArrNum="$OPTARG";;
        p) AddPrm="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`
PartFil=`readlink -f $PartFil`

FamNam=`cut -f 1 $PartFil | head -n $ArrNum | tail -n 1`
Proband=`cut -f 2 $PartFil | head -n $ArrNum | tail -n 1`
Parent=`cut -f 3 $PartFil | head -n $ArrNum | tail -n 1`
Extras=`cut -f 4 $PartFil | head -n $ArrNum | tail -n 1`


echo $FamNam
if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.AR --alt $Proband --het $Parent"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.AR.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.AR.tsv &
fi
#Autosomal Dominant unaffected Parent
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.unaffAD  --het $Proband --ref $Parent"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.unaffAD.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.unaffAD.tsv &
fi
#Autosomal Dominant affected Parent
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.affAD  --het $Proband,$Parent"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.affAD.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.affAD.tsv &
fi
#compound heterozygous
echo "Compund heterozygous.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.tempheppat  --het $Proband,$Parent -P  -f 0.03"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.temphepmat  --het $Proband --ref $Parent -P  -f 0.03"
if [[ -n $Extras ]]; then CMD=$CMD" --unfl $Extras"; fi
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$FamNam.ParentChild.temphepmat.tsv")
pathet <- read.delim("$FamNam.ParentChild.tempheppat.tsv")

mathet <- mathet[,match(colnames(pathet), colnames(mathet))]

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$FamNam.ParentChild.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
cat $FamNam.ParentChild.tempheppat.log $FamNam.ParentChild.temphepmat.log > $FamNam.ParentChild.compound_heterozygous.log
LEN=`cat $FamNam.ParentChild.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.compound_heterozygous.tsv &
fi
rm -rf *temp*
