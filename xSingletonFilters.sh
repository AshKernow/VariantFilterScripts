#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterSing

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n:p: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) SingFil="$OPTARG";;
        p) AddPrm="$OPTARG";;
        n) DirPre="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
SingFil=`readlink -f $SingFil`

Proband=`cut -f 2 $SingFil | head -n $SGE_TASK_ID | tail -n 1`
FamNam=`cut -f 1 $SingFil | head -n $SGE_TASK_ID | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $Proband.AR --alt $Proband"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $Proband.AR.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $Proband.AR.tsv
fi

#Autosomal Dominant
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $Proband.AD  --het $Proband"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $Proband.AD.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $Proband.AD.tsv
fi
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

het <- read.delim("$Proband.AD.tsv")

gens <- unique(het[duplicated(het[,"Gene"]),"Gene"])

comphet <- het[het[,"Gene"]%in%gens,]

write.table(comphet, "$Proband.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
LEN=`cat $Proband.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $Proband.compound_heterozygous.tsv
fi
