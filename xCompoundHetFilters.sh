#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterCHet

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n:p: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) FamFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        p) AddPrm="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
FamFil=`readlink -f $FamFil`

FamNam=`cut -f 1 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
ModNam=`cut -f 2 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
PatNam=`cut -f 3 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
MatNam=`cut -f 4 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
FilPrm=`cut -f 5 $FamFil | head -n $SGE_TASK_ID | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam.Fam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

PatPrm=`echo $FilPrm | sed s/"--het "/"--het $PatNam,"/`
if [[ "$PatPrm" == *--ref* ]]; then
    PatPrm=`echo $PatPrm | sed s/"--ref "/"--ref $MatNam,"/`
else
    PatPrm=$PatPrm" --ref "$MatNam
fi
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.tempheppat $PatPrm"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
MatPrm=`echo $FilPrm | sed s/"--het "/"--het $MatNam,"/`
if [[ "$MatPrm" == *--ref* ]]; then
    MatPrm=`echo $MatPrm | sed s/"--ref "/"--ref $MatNam,"/`
else
    MatPrm=$MatPrm" --ref "$PatNam
fi
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.temphepmat $MatPrm"
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
rm -rf *temphep*

