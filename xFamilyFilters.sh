#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterFam

ArrNum=1

usage="xFamilyFilters.sh -v <vcf file> -t <Filter Table> -n <Output Directory> -a <The line of the Trio Table to analyse> -p <Additional Parameters> -H <this message>"

#get arguments
while getopts v:t:n:a:p:H opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) FamFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        a) ArrNum="$OPTARG";;
        p) AddPrm="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`
FamFil=`readlink -f $FamFil`


FamNam=`cut -f 1 $FamFil | head -n $ArrNum | tail -n 1`
ModNam=`cut -f 2 $FamFil | head -n $ArrNum | tail -n 1`
FilPrm=`cut -f 3 $FamFil | head -n $ArrNum | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#check for X-linked and add flag if necessary
if [[ $ModNam == *X* ]]; then FilPrm=$FilPrm" -X"; fi
#check for de novo and add flag if necessary
if [[ $ModNam == *ovo* ]]; then FilPrm=$FilPrm" -D"; fi
#check for Individual compound heterozygous and flag if necessary
if [[ "$ModNam" == *IcHet* ]]; then FilPrm=$FilPrm" -P"; fi

echo "Filtering.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.$ModNam $FilPrm"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
if [[ "$ModNam" == *IcHet* ]]; then CMD=$CMD" -f 0.03"; fi
echo $CMD
eval $CMD

#check for Individual compound heterozygous and 
if [[ "$ModNam" == *IcHet* ]]; then
    TsvFil=$FamNam.$ModNam.tsv
    LogFil=$FamNam.$ModNam.log
    echo "----------------------------------------------------------------------" >> $LogFil
    echo "Filtering for compound heterozygous..." >> $LogFil
    R --vanilla <<RSCRIPT
        dat <- read.delim("$TsvFil")
        dat <- dat[gsub(",.*", "", dat[,6])%in%gsub(",.*", "", dat[duplicated(dat[,6]),6]),]
        dat <- dat[dat[,6]!=".",]
        write.table(dat, "$TsvFil", col.names=T, row.names=F, sep="\t", quote=F)
RSCRIPT
    LEN=`cat $TsvFil | wc -l`
    LEN=$(( LEN - 1 ))
    GEN=`tail -n +2 $TsvFil | cut -f 6 | uniq | wc -l`
    echo "After filtering $LEN variants remaining in $GEN genes" >> $LogFil
fi
LEN=`cat $FamNam.$ModNam.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.$ModNam.tsv &
fi
