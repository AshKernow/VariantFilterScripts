#!/bin/bash
#$ -l mem=4G,time=6:: -cwd -S /bin/bash -N SplitbyFam

#set default arguments
usage="
ExmVCFtoFamilies.sh -i <InputFamilyTable> -v <VcfFile> -o <OutputPrefixName>

     -i (required) - Table with famly list <FAMILY_ID>\t<FAMILY_SAMPLE_LIST>; family sample list is comma separated
     -v (required) - VCFfile to be split
     -o (optional) - A prefix to add to the output file names
     -K (flag) - Keep all variants; default is to only keep variants with an alternate allele in the new vcf
     -T (flag) - Convert to tabular format
     -H (flag) - echo this message and exit
"

#get arguments
AltAll="true"
ToTable="false"
while getopts i:v:o:KTH opt; do
    case "$opt" in
        i) InpTab="$OPTARG";;
        v) VcfFil="$OPTARG";;
        o) OutPre="$OPTARG";;
        K) AltAll="false";;
        T) ToTable="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
InpTab=`readlink -f $InpTab`
VcfFil=`readlink -f $VcfFil`
if [[ ! -e "$InpTab" ]] | [[ ! -e "$VcfFil" ]]; then 
    echo "Missing/Incorrect required arguments"
    echo "Family Table: $InpTab; VCF file: $VcfFil"
    echo "$usage"
    exit $1
fi
RowNum=$SGE_TASK_ID
if [[ "$RowNum" == "undefined" ]]; then RowNum=1; fi
FamNam=`tail -n +$RowNum $InpTab | head -n 1 | cut -f 1`
FamSam=`tail -n +$RowNum $InpTab | head -n 1 | cut -f 2 `
FamPrm=`echo $FamSam | sed s/\"//g | sed s/\ //g | sed s/[*mf]$// | sed s/[*mf],/,/g  | sed s/,/\ --indv\ /g`
FamPrm=" --indv "$FamPrm

echo "Input table: "$InpTab"; Row: "$RowNum
echo "Family: "$FamSam" Samples:"$FamSam
echo 

OutVcf=`basename $VcfFil`
OutVcf=${OutVcf/.vcf/}"_"$OutPre$FamNam
TmpVCF=TEMP$OutVcf

CMD="vcftools --vcf $VcfFil $FamPrm --recode --recode-INFO-all --out $TmpVCF"
echo $CMD
eval $CMD

if [[ $AltAll == "true" ]]; then 
    CMD="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/ExmFilt.3.AlternateAlleleFilter.py -v $TmpVCF.recode.vcf -o $OutVcf"
    echo
    echo $CMD
    eval $CMD
else
    mv $TmpVCF.recode.vcf $OutVcf.vcf
fi
rm -f $TmpVCF*

if [[ "$ToTable" == "true" ]]; then
    CMD="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/ExmFilt.4.VCFtoTable.py -v $OutVcf.vcf -o $OutVcf"
    echo
    echo $CMD
    eval $CMD
fi
