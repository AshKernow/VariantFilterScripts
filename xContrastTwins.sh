#!/bin/bash
#$ -cwd -l mem=1G,time=:30: -N FilterTwin

#get arguments
while getopts v:t:n:p:f:a:b:u: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) TwnTab="$OPTARG";;
        n) DirPre="$OPTARG";;
        p) AddPrm="$OPTARG";;
        f) FamNam="$OPTARG";;
        a) Twin1="$OPTARG";;
        b) Twin2="$OPTARG";;
        u) Unfiltered="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`

if [[ $TwnTab ]]; then
    FamNam=`cut -f 1 $TwnTab | head -n $SGE_TASK_ID | tail -n 1`
    Twin1=`cut -f 2 $TwnTab | head -n $SGE_TASK_ID | tail -n 1`
    Twin2=`cut -f 3 $TwnTab | head -n $SGE_TASK_ID | tail -n 1`
    Unfiltered=`cut -f 4 $TwnTab | head -n $SGE_TASK_ID | tail -n 1`
fi

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi
OutFil=$FamNam.ContrastTwins

DirNam=$FamNam
if [[ $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

runCustomFilter(){
if [[ $Unfiltered ]]; then TwinComp=$TwinComp" --unfl $Unfiltered"; fi
if [[ $AddPrm ]]; then TwinComp=$TwinComp" $AddPrm"; fi
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil $TwinComp"
echo $CMD
eval $CMD
}

TwinComp=" -o $OutFil.t1.temp --het $Twin1 --ref $Twin2"
runCustomFilter

TwinComp=" -o $OutFil.t2.temp --het $Twin1 --alt $Twin2"
runCustomFilter

TwinComp=" -o $OutFil.t3.temp --ref $Twin1 --het $Twin2"
runCustomFilter

TwinComp=" -o $OutFil.t4.temp --ref $Twin1 --alt $Twin2"
runCustomFilter

TwinComp=" -o $OutFil.t5.temp --alt $Twin1 --ref $Twin2"
runCustomFilter

TwinComp=" -o $OutFil.t6.temp --alt $Twin1 --het $Twin2"
runCustomFilter

R --vanilla <<RSCRIPT
options(stringsAsFactors=F)
dat <- read.delim("$OutFil.t1.temp.tsv")
for(i in 2:6){
    temp <- read.delim(paste("$OutFil.t", i, ".temp.tsv", sep=""))
    temp <- temp[,match(colnames(temp), colnames(dat))]
    dat <- rbind(dat, temp)
    }
write.table(dat, "$OutFil.tsv", row.names=F, quote=F, sep="\t")
RSCRIPT

LEN=`cat $OutFil.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $OutFil.tsv
fi


echo "Six contrasts performed on paired samples. Logs for each contrast:" > $OutFil.filtered.log
echo >> $OutFil.filtered.log
echo "-------------------------------------------------------------------------------------------------------------" >> $OutFil.filtered.log
echo >> $OutFil.filtered.log

for i in $OutFil.t*.temp.log; do 
 grep -vE  "OutputName|Not" $i >> $OutFil.log
 echo >> $OutFil.filtered.log
 echo "-------------------------------------------------------------------------------------------------------------" >> $OutFil.filtered.log
 echo >> $OutFil.filtered.log
done

rm -f $OutFil.*temp*
