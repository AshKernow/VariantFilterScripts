#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N AnnTSV

# The purpose of this script to annotate a tab separated files containing lists of target variants. The file firsts gets ClinVar annotations from Annovar. The R script adds annotations from ClinVar at the variants level, and then annotations from a custom annotation table at the gene level.

usage="-i <InputFile> -o <OutputName - optional> -H <This message>"

#get arguments
while getopts i:o:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        o) OutFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#if not output name has been specified derive it from the input file name
if [[ -z $OutFil ]]; then 
    OutFil=${InpFil/tsv/annot.tsv}
else
    OutFil=$OutFil.annot.tsv
fi

# location of annovar database adn the custom annotation table
AnnDir=/home/local/ARCS/ads2202/resources/annovar
AnnTab=/home/local/ARCS/ads2202/resources/VariantAnnotationTables/Final_Gene_Annotation_Table.txt

echo $InpFil
echo $OutFil

#check input has variants to annotate
LinNum=`cat $InpFil | wc -l`
if [[ $LinNum -lt 2 ]];then
    echo "There are no variants to annotate"
    exit
fi


TempInp=$InpFil.tempann #name for temporary files

#cut the first 8 column of the tsv to simulate vcf format
tail -n +2 $InpFil | cut -f 1-8 > $TempInp

#convert to annovar format - annovar adjusts the way indels are represented
CMD1="convert2annovar.pl -includeinfo -format vcf4 $TempInp > $TempInp.2; mv -f $TempInp.2 $TempInp"

#annotate the new file
CMD2="table_annovar.pl $TempInp $AnnDir/humandb/ -buildver hg19 -remove -protocol clinvar_20150330 -operation f -nastring . -otherinfo -outfile $TempInp"

#Run commands
echo "Start Annovar at `date`"
echo $CMD1
eval $CMD1
echo $CMD2
eval $CMD2


##annotate using R
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

#get data
dat.all <- read.delim("$InpFil", colClasses="character")
dat <- dat.all[,1:(grep("AlternateAlleles", colnames(dat.all)))]

#get the population specific frequencies
popfrqcol <- which(colnames(dat.all)%in%c("ESP.aa", "ESP.ea", "1KG.eur", "1KG.amr", "1KG.eas", "1KG.afr", "1KG.sas", "ExAC.afr", "ExAC.amr", "ExAC.eas", "ExAC.fin", "ExAC.nfe", "ExAC.oth", "ExAC.sas"))
POPFREQ <- dat.all[,popfrqcol]

# Add a AAF column for Wen-I
ADcols <- grep("AD$", colnames(dat))
GTcols <- grep("GT$", colnames(dat))
for(j in 1:length(ADcols)){
    dat <- cbind(dat, 0)
    ADcol <- ADcols[j]
    GTcol <- GTcols[j]
    nam <- gsub("AD", "AAF", colnames(dat)[ADcol])
    colnames(dat)[ncol(dat)] <- nam
    for(i in 1:nrow(dat)){
        AAs <- as.numeric(strsplit(dat[i,GTcol], "/")[[1]])
        AAs <- unique(AAs[AAs!=0]+1)
        ADs <- as.numeric(strsplit(dat[i,ADcol], ",")[[1]])
        tot <- sum(ADs)
        AAFs <- round(ADs[AAs]/tot,5)
        if(length(AAFs)==0) AAFs <- 0
        dat[i,nam] <- paste(AAFs, collapse=",")
    }
}
#add leading space to genotypes to deal with Excel date format issue and read depths to stop auto correction to numbers
for(i in c(GTcols, ADcols)){
    dat[,i] <- paste(" ", dat[,i], sep="")
    dat[,i] <- gsub(",", ", ", dat[,i])
}
#get annotation table and extract relevant lines
annot <- read.delim("$AnnTab")
annmat <- match(gsub(",.*", "", dat[,"Gene"]), annot[,"GENE"])
annot <- annot[annmat,]
if(any(is.na(annot[,1]))) { annot[is.na(annot[,1]),] <- "." } #replace missing data (NA) with "."

#get clinvar and match by variatns (CHR_POS_REF_ALT)
CLINVAR <- read.table("$TempInp.hg19_multianno.txt", skip=1)
  #need to match to specific variant, use the original variants that are at the end of the annovar output
clinmat <- match(paste(dat[,1], dat[,2], dat[,4], dat[,5]), paste(CLINVAR[,7], CLINVAR[,8], CLINVAR[,10], CLINVAR[,11]))
CLINVAR <- CLINVAR[clinmat,6]
if(any(is.na(CLINVAR))){CLINVAR[is.na(CLINVAR)] <- "."} #replace missing data (NA) with "."

#Create TOLERANCE.specific column with variant types
TOLERANCE.specific <- rep(".", nrow(dat))
whi <- which(dat[,"VariantClass"]=="synonymousSNV")
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_SYNONYMOUS"]}
whi <- which(dat[,"VariantClass"]=="nonsynonymousSNV")
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_MISSENSE"]}
whi <- which(dat[,"VariantClass"]%in%c("stopgain", "stoploss"))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_NONSENSE"]}
whi <- which(grepl("splicing", dat[,"VariantFunction"]))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_SPLICE"]}
whi <- which(grepl("^frameshift", dat[,"VariantClass"]))
if(length(whi)>0){TOLERANCE.specific[whi] <- annot[whi,"TOLERANCE_FRAMESHIFT"]}

#insert TOLERANCE.specific column into annot and remove individual columns for each variant type
col1 <- grep("TOLERANCE_ALL_DALY", colnames(annot))
col2 <- grep("TOLERANCE_FRAMESHIFT", colnames(annot))+1
ANNOT <- cbind(annot[,2:col1], TOLERANCE.specific, annot[,col2:ncol(annot)])

#combine annotations
out <- cbind(dat, CLINVAR, ANNOT, POPFREQ)


#sort the file so that most likely candidates are at the top...
predord <- match(out[,"PredictionSummary"], c("Med", "High"))        #...sort high likelihood deleterious missense to the top
splicord <- grepl("splicing", out[,"VariantFunction"])               #...then splicing to the nonsense to the top of those
nonseord <- grepl("stop", out[,"VariantClass"])                      #...then stopgain/stoploss
indelfsord <- grepl("^frameshift", out[,"VariantClass"])             #...then frameshift indels
indelnford <- grepl("nonframeshift", out[,"VariantClass"])           #...then non-frameshift indels
qualord <- match(out[,"VariantCallQuality"], c("Low","Med", "High")) #...then high quality variants
cadd <- out[,"CADDscore"]                                            #...final sort by CADD and then GERP
gerp <- out[,"GERP.."]                                               #...and then GERP
ord <- order(predord, splicord, nonseord, indelfsord, indelnford, qualord, cadd, gerp, decreasing=T)
out <- out[ord,]

#Truncate fields that are too long for a single Excel cell
for(i in 1:ncol(out)){
  cnt <- which(nchar(out[,i])>32760)
  if(length(cnt)>0){
    out[cnt,i] <- paste(substr(out[cnt,i], 1, 32700), ":::---TRUNCATED_FOR_EXCEL")
  }
}

## Adjust column headers
colnames(out) <- gsub("GERP..", "GERP++", colnames(out))
colnames(out) <- gsub("prediction", "", colnames(out))
colnames(out) <- gsub("^MA$", "Mutation Assessor", colnames(out))
colnames(out) <- gsub("^MT$", "Mutation Taster", colnames(out))
colnames(out) <- gsub("^PP2$", "PolyPhen2", colnames(out))
colnames(out) <- gsub("^CADDscore$", "CADD score", colnames(out))

#OUTPUT
write.table(out, "$OutFil", sep="\t", col.names=T, row.names=F, quote=F)
RSCRIPT

if [[ $? -gt 0 ]]; then
    echo "Annotation Failed"
    exit
fi

# remove temporary files
rm -f $InpFil.tempann*
echo "Annotation complete"
