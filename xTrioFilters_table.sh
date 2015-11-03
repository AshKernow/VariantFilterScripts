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

TempFiltScrDir="/home/local/ARCS/ads2202/scripts/TESTINGFILTER"
FiltScrDir="/home/local/ARCS/ads2202/scripts/Filtering_scripts"

VcfFil=`readlink -f $VcfFil`
TrioFil=`readlink -f $TrioFil`


FamNam=`cut -f 1 $TrioFil | head -n $ArrNum | tail -n 1`
Proband=`cut -f 2 $TrioFil | head -n $ArrNum | tail -n 1`
Father=`cut -f 3 $TrioFil | head -n $ArrNum | tail -n 1`
Mother=`cut -f 4 $TrioFil | head -n $ArrNum | tail -n 1`
Unfiltered=`cut -f 5 $TrioFil | head -n $ArrNum | tail -n 1`


echo $FamNam
if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#AR, AD
FilterTable=$FamNam.Trio.MainFilterTable.tab
echo -e "#FAMILY\tMODEL\tSamples_0/0\tSamples_0/1\tSamples_1/1\tSamples_not_0/0\tSamples_1/1\tUnfiltered_Samples" > $FilterTable
echo -e "$FamNam\tAR\t\t$Father,$Mother\t$Proband\t\t\t$Unfiltered" >> $FilterTable
echo -e "$FamNam\tAD-paternal\t$Mother\t$Proband,$Father\t\t\t\t$Unfiltered" >> $FilterTable
echo -e "$FamNam\tAD-maternal\t$Father\t$Proband,$Mother\t\t\t\t$Unfiltered" >> $FilterTable
CMD="$TempFiltScrDir/ExmFilt.CustomGenotype_with_Table.py -v $VcfFil -t $FilterTable -o Trio"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
for Model in AR AD-paternal AD-maternal; do
    LEN=`cat $FamNam.Trio.$Model.tsv | wc -l`
    if [[ $LEN -gt 1 ]]; then
        nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.$Model.tsv &
    fi
done

#de novo
FilterTable=$FamNam.Trio.DenovoFilterTable.tab
echo -e "#FAMILY\tMODEL\tSamples_0/0\tSamples_0/1\tSamples_1/1\tSamples_not_0/0\tSamples_1/1\tUnfiltered_Samples" > $FilterTable
echo -e "$FamNam\tdenovo\t$Father,$Mother\t$Proband\t\t\t\t$Unfiltered" >> $FilterTable
CMD="$TempFiltScrDir/ExmFilt.CustomGenotype_with_Table.py -v $VcfFil -t $FilterTable -o Trio"
if [[ "$BadDeN" == "false" ]]; then CMD=$CMD" -D"; fi #otherwise de novos will be run with default filters
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.denovo.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.denovo.tsv &
fi

# X-linked
FilterTable=$FamNam.Trio.X-linkedFilterTable.tab
echo -e "#FAMILY\tMODEL\tSamples_0/0\tSamples_0/1\tSamples_1/1\tSamples_not_0/0\tSamples_1/1\tUnfiltered_Samples" > $FilterTable
echo -e "$FamNam\tX-linked\t$Father\t$Proband,$Mother\t\t\t\t$Unfiltered" >> $FilterTable
CMD="$TempFiltScrDir/ExmFilt.CustomGenotype_with_Table.py -v $VcfFil -t $FilterTable -o Trio -X"
if [[ ! -z $AddPrm ]]; then CMD=$CMD" $AddPrm"; fi
echo $CMD
eval $CMD
LEN=`cat $FamNam.Trio.X-linked.tsv | wc -l`
if [[ $LEN -gt 1 ]]; then
    nohup $FiltScrDir/xAnnotateVariantTSV.sh -i $FamNam.Trio.X-linked.tsv &
fi

# Compound het
FilterTable=$FamNam.Trio.tempcHetFilterTable.tab
echo -e "#FAMILY\tMODEL\tSamples_0/0\tSamples_0/1\tSamples_1/1\tSamples_not_0/0\tSamples_1/1\tUnfiltered_Samples" > $FilterTable
echo -e "$FamNam\ttempheppat\t$Mother\t$Proband,$Father\t\t\t\t$Unfiltered" >> $FilterTable
echo -e "$FamNam\ttemphepmat\t$Father\t$Proband,$Mother\t\t\t\t$Unfiltered" >> $FilterTable
CMD="$TempFiltScrDir/ExmFilt.CustomGenotype_with_Table.py -v $VcfFil -t $FilterTable -o Trio -P -f 0.03"
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


#rm -rf $FamNam.Trio.temp*
