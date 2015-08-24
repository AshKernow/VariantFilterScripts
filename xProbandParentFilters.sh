#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterPC

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n:p: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) PartFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        p) AddPrm="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
PartFil=`readlink -f $PartFil`

FamNam=`cut -f 1 $PartFil | head -n $SGE_TASK_ID | tail -n 1`
Proband=`cut -f 2 $PartFil | head -n $SGE_TASK_ID | tail -n 1`
Parent=`cut -f 3 $PartFil | head -n $SGE_TASK_ID | tail -n 1`


echo $FamNam
if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.AR --alt $Proband --het $Parent"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.AR.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.AR.tsv
fi
#Autosomal Dominant unaffected Parent
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.unaffAD  --het $Proband --ref $Parent"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.unaffAD.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.unaffAD.tsv
fi
#Autosomal Dominant affected Parent
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.affAD  --het $Proband,$Parent"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.ParentChild.affAD.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.affAD.tsv
fi
#compound heterozygous
echo "Compund heterozygous.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.tempheppat  --het $Proband,$Parent -P  -f 0.03"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.ParentChild.temphepmat  --het $Proband --ref $Parent -P  -f 0.03"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$FamNam.ParentChild.temphepmat.tsv")
pathet <- read.delim("$FamNam.ParentChild.tempheppat.tsv")

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$FamNam.ParentChild.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
cat $FamNam.ParentChild.tempheppat.log $FamNam.ParentChild.temphepmat.log > $FamNam.ParentChild.compound_heterozygous.log
LEN=`cat $FamNam.ParentChild.compound_heterozygous.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    qsub $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.ParentChild.compound_heterozygous.tsv
fi
rm -rf *temp*
