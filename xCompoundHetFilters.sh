#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterCHet


## use "0" for missing parent
ArrNum=1
#get arguments
while getopts v:t:n:a:p: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) FamFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        a) ArrNum="$OPTARG";;
        p) AddPrm="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

echo "Starting"

FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`
FamFil=`readlink -f $FamFil`

FamNam=`cut -f 1 $FamFil | head -n $ArrNum | tail -n 1`
ModNam=`cut -f 2 $FamFil | head -n $ArrNum | tail -n 1`
PatNam=`cut -f 3 $FamFil | head -n $ArrNum | tail -n 1`
MatNam=`cut -f 4 $FamFil | head -n $ArrNum | tail -n 1`
FilPrm=`cut -f 5 $FamFil | head -n $ArrNum | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

LogFil=$FamNam.compound_heterozygous.log

echo "Start :`date`"
echo "VCF file: $VcfFil"
echo "Filter file: $FamFil"
if [[ ! -e $VcfFil ]] | [[ ! -e $FamFil ]]; then "Missing/Incorrect required arguments"; exit; fi




#get paternally interited allele
PatPrm=$FilPrm
if [[ "$PatNam" != "0" ]]; then 
    PatPrm=`echo $PatPrm | sed s/"--het "/"--het $PatNam,"/`
fi
if [[ "$MatNam" != "0" ]]; then 
    if [[ "$PatPrm" == *--ref* ]]; then
        PatPrm=`echo $PatPrm | sed s/"--ref "/"--ref $MatNam,"/`
    else
        PatPrm=$PatPrm" --ref "$MatNam
    fi
fi
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.tempheppat $PatPrm -P -f 0.03"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
echo `date`
eval $CMD


MatPrm=$FilPrm
if [[ "$MatNam" != "0" ]]; then 
    MatPrm=`echo $MatPrm | sed s/"--het "/"--het $MatNam,"/`
fi
if [[ "$PatNam" != "0" ]]; then 
    if [[ "$MatPrm" == *--ref* ]]; then
        MatPrm=`echo $MatPrm | sed s/"--ref "/"--ref $PatNam,"/`
    else
        MatPrm=$MatPrm" --ref "$PatNam
    fi
fi
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.temphepmat $MatPrm -P -f 0.03"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
echo `date`
eval $CMD

echo "Run R script to find compound het"
echo `date`

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

echo "Paternal Filtering" > $LogFil
echo "---------------------------------------------------" >> $LogFil
cat $FamNam.tempheppat.log >> $LogFil
echo "" >> $LogFil
echo "===================================================" >> $LogFil
echo "" >> $LogFil

echo "Maternal Filtering" >> $LogFil
echo "---------------------------------------------------" >> $LogFil
cat $FamNam.temphepmat.log >> $LogFil
echo "" >> $LogFil
echo "===================================================" >> $LogFil
echo "" >> $LogFil


rm -rf *temphep*
LEN=`cat $FamNam.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    CMD="nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.compound_heterozygous.tsv &"
    echo $CMD >> $LogFil
    eval $CMD
fi

GENLEN=`cut -f 6 $FamNam.compound_heterozygous.tsv | tail -n +2 | sort | uniq | wc -l`


echo "Final Number of variants: $LEN in $GENLEN genes" >> $LogFil
echo "" >> $LogFil
echo "===================================================" >> $LogFil

echo "Finished `date`"