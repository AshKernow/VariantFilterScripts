#!/usr/bin/env python
#$ -N CustFilt -cwd

################################################################################################################
##### USAGE                                                                                                  ###
################################################################################################################
#### The purpose of this script is to tabulare a VCF file, one line for each variant allele, all of the interesting INFO fields in columns
####    -v/--vcf    <required>    The script requires an input vcf
####    -o/--out    <required>    The user should specify a base name for the output files.
####    -p/--pos ...
################################################################################################################
################################################################################################################

################################################################################################################
##### GET PARAMETERS                                                                                         ###
################################################################################################################
import os
from optparse import OptionParser
parser = OptionParser()

## Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="Input VCF file", metavar="VCFfile")
parser.add_option("-o", "--out", dest="OutputFileName",help="Base of names for output files", metavar="OutputFileName")

(options, args) = parser.parse_args()


################################################################################################################
##### ASSIGN FILTER VARIABLES                                                                                ###
################################################################################################################

## Variant Quality Filters
BadSnpFilters=['QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP']
BadInDFilters=['LowQD_Indel']
MidSnpFilters=['FS_Mid_SNP','QD_Mid_SNP','VQSRTrancheSNP99.00to99.90','VQSRTrancheSNP99.90to100.00']
MidInDFilters=['FSBias_Indel','RPBias_Indel','VQSRTrancheINDEL99.00to99.90','VQSRTrancheINDEL99.90to100.00']

## Variant classes
MissenseClass=['nonsynonymousSNV','unknown']
NonFrameshiftClass=['nonframeshiftdeletion','nonframeshiftinsertion','nonframeshiftsubstitution']
FrameshiftClass=['frameshiftdeletion','frameshiftinsertion','frameshiftsubstitution']
InDelClass=NonFrameshiftClass+FrameshiftClass
NonsenseClass=['stopgain','stoploss']
SplicingClass=['splicing','exonic,splicing']
################################################################################################################
##### OPEN INPUT FILES AND OUTPUT FILES                                                                      ###
################################################################################################################

##Assign input and output files
VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
TabOutputFilename=BaseName+'.tsv'
Output=open(TabOutputFilename,'w')
LogOutputFilename=BaseName+'.log'
Outlog=open(LogOutputFilename,'w')

##write log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Tabulate VCF: "+TimeNow+"\n")
Outlog.write("VCF: "+str(options.VCFfile)+"\n")
Outlog.write("OutputName: "+str(options.OutputFileName)+"\n")
Outlog.write("\n")

################################################################################################################
##### START FILTERING                                                                                        ###
################################################################################################################

OrigCount=0
for line in VCF:
    if  line.startswith('#CHROM'):
        ################################################################################################################
        ##### INITIALISE OUTPUT FILES                                                                                ###
        ################################################################################################################
        linelist=line.split("\t", 9)
        SampleList=linelist[9].strip()
        AllSamplesList=SampleList.split("\t")
        ##Start output tables
        headerlist=['Chromosome','Position','ID','REF','ALT', 'AlleleNum', 'Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.ExAC','AlleleFrequency.1KG','AlleleFrequency.ESP','AlleleFrequency.VCF','SIFTprediction','PP2prediction','MAprediction','MTprediction','GERP++','CADDscore','MetaSVMprediction','SegmentalDuplication','PredictionSummary','VariantCallQuality','QUAL','MQ0','DP','FILTER']+[ i+" GT" for i in AllSamplesList]+[ i+" AD" for i in AllSamplesList]+[ i+" DP" for i in AllSamplesList]+[ i+" GQ" for i in AllSamplesList]+['AlternateAlleles', 'INFO']
        Output.write("\t".join(headerlist)+"\n")
    if '#' not in line:
        ################################################################################################################
        ##### PARSE LINE AND GET VARIANT INFO                                                                        ###
        ################################################################################################################
        line=line.strip()
        linelist=line.split("\t")
        OrigCount=OrigCount+1
        
        ##Get Alternate Alleles
        AltAllsStr=linelist[4]
        AltAlls=AltAllsStr.split(",")
        AltNum=len(AltAlls)
        
        QUAL=float(linelist[5])
        VariantFilter=linelist[6]
        FORMAT=linelist[8]
        FORMAT=FORMAT.split(":")
        
        ## Get variant data from INFO Field
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        DPnumber=float(INFOdict.get('DP','.'))
        MQ0number=float(INFOdict.get('MQ0','.'))
        GeneName=INFOdict.get('GeneName','.')
        VariantFunction=INFOdict.get('VarFunc','none')
        SegDup=INFOdict.get('SegDup','none')
        
        KGFreqList=str(INFOdict.get('1KGfreq','.'))
        KGFreqList=KGFreqList.split(',')
        ESPFreqList=str(INFOdict.get('ESPfreq','.'))
        ESPFreqList=ESPFreqList.split(',')
        ExACFreqList=str(INFOdict.get('ExACfreq','.'))
        ExACFreqList=ExACFreqList.split(",")
        VCFFreqList=str(INFOdict.get('AF',0))
        VCFFreqList=VCFFreqList.split(",")
        
        VariantClassList=INFOdict.get('VarClass','none')
        VariantClassList=VariantClassList.split(',')
        AAchangeList=INFOdict.get('AAChange','.')
        AAchangeList=AAchangeList.split(',')
        SIFTpredictionList=INFOdict.get('SIFTprd','.')
        SIFTpredictionList=SIFTpredictionList.split(',')
        PP2predictionList=INFOdict.get('PP2.hvar.prd','.')
        PP2predictionList=PP2predictionList.split(',')
        MApredictionList=INFOdict.get('MutAprd','.')
        MApredictionList=MApredictionList.split(',')
        MTpredictionList=INFOdict.get('MutTprd','.')
        MTpredictionList=MTpredictionList.split(',')
        GERPscoreList=str(INFOdict.get('GERP','.'))
        GERPscoreList=GERPscoreList.split(',')
        CADDscoreList=str(INFOdict.get('CADDphred','.'))
        CADDscoreList=CADDscoreList.split(',')
        MetaSVMpredictionList=INFOdict.get('MetaSVMprd','.')
        MetaSVMpredictionList=MetaSVMpredictionList.split(',')
        
        
        
        ################################################################################################################
        ##### CHECK TO SEE IF GENOTYPES PASS, ITERATE ACROSS MULTIPLE ALT ALLELES IF NECESSARY                      ###
        ################################################################################################################
        
        AltRng=range(0, AltNum)
        for altnum in AltRng:
            
            ################################################################################################################
            ##### VARIANT SPECIFIC CLASS CHECK, FREQUENCY AND VARIANT CALLING QUALITY                                    ###
            ################################################################################################################
            
            ##check variant class
            cltnum=min(len(VariantClassList)-1, altnum)
            VariantClass=str(VariantClassList[cltnum])
            
            ## Check if KG passes threshold
            cltnum=min(len(KGFreqList)-1, altnum)
            KGFreq=str(KGFreqList[cltnum])
            
            ## Check if ESP passes threshold
            cltnum=min(len(ESPFreqList)-1, altnum)
            ESPFreq=str(ESPFreqList[cltnum])
                
            ## Check if ExAC passes threshold
            cltnum=min(len(ExACFreqList)-1, altnum)
            ExACFreq=str(ExACFreqList[cltnum])
            
            ## Check if VCF passes threshold
            cltnum=min(len(VCFFreqList)-1, altnum)
            VCFFreq=str(VCFFreqList[cltnum])
            
            ## Check and annotate pathogenicity (discard missense predicted to be benign - need to adjust this in splicining regions)
            
            cltnum=min(len(SIFTpredictionList)-1, altnum)
            SIFTprediction=SIFTpredictionList[cltnum]
            cltnum=min(len(PP2predictionList)-1, altnum)
            PP2prediction=PP2predictionList[cltnum]
            cltnum=min(len(MetaSVMpredictionList)-1, altnum)
            MetaSVMprediction=MetaSVMpredictionList[cltnum]
            cltnum=min(len(CADDscoreList)-1, altnum)
            CADDscore=str(CADDscoreList[cltnum])
            if CADDscore == ".":
                CADDscoretest=float(0)
            else:
                CADDscoretest=float(CADDscore)
            
            #Set Patho level
            PathoLevel="Low"
            if VariantClass in MissenseClass and SIFTprediction=="." and PP2prediction=="." and CADDscore=="." and MetaSVMprediction==".":
                PathoLevel="Med"
            if VariantClass in MissenseClass and (SIFTprediction=="D" or PP2prediction=="D" or PP2prediction=="P" or CADDscoretest>=15):
                PathoLevel="Med"
            if VariantClass in MissenseClass and (SIFTprediction=="D" and PP2prediction=="D"):
                PathoLevel="High"
            if VariantClass in MissenseClass and (CADDscoretest>=25):
                PathoLevel="High"
            if VariantClass in MissenseClass and MetaSVMprediction=="D":
                PathoLevel="High"
            if VariantFunction in SplicingClass:
                PathoLevel="High"
            if VariantClass in InDelClass or VariantClass in NonsenseClass:
                PathoLevel="High"
            
            ##Set Filter Quality summary
            FILTER='High'
            if VariantClass in InDelClass and any( str(i) in VariantFilter for i in MidInDFilters):
                FILTER='Medium'
            if VariantClass not in InDelClass and any( str(i) in VariantFilter for i in MidSnpFilters):
                FILTER='Medium'
            if VariantClass in InDelClass and any( str(i) in VariantFilter for i in BadInDFilters):
                FILTER='Low'
            if VariantClass not in InDelClass and any( str(i) in VariantFilter for i in BadSnpFilters):
                FILTER='Low'
            
            ################################################################################################################
            ##### OUTPUTS                                                                                                ###
            ################################################################################################################
            
            ## If all pass then output line
            
            cltnum=min(len(AAchangeList)-1, altnum)
            AAchange=AAchangeList[cltnum]
            cltnum=min(len(MApredictionList)-1, altnum)
            MAprediction=MApredictionList[cltnum]
            cltnum=min(len(MTpredictionList)-1, altnum)
            MTprediction=MTpredictionList[cltnum]
            ##output
            ALT=str(AltAlls[altnum])
            cltnum=min(len(GERPscoreList)-1, altnum)
            GERPscore=str(GERPscoreList[cltnum])
            SampleStrings=linelist[9:]
            SampleStrings=[ i.split(':') for i in SampleStrings ]
            GTind=FORMAT.index('GT')
            AllSampleGT=[ "." for i in range(0,len(SampleStrings))]
            for i in range(0,len(SampleStrings)):
                AllSampleGT[i]=SampleStrings[i][GTind]
            
            if "AD" in FORMAT:
            ## Define Allele Count
                ADind=FORMAT.index('AD')
                AllSampleAD=[ "." for i in range(0,len(SampleStrings))]
                for i in range(0,len(SampleStrings)):
                    if AllSampleGT[i]!="./.":
                        AllSampleAD[i]=SampleStrings[i][ADind]
            
            if "GQ" in FORMAT:
            ## Define Allele Count
                GQind=FORMAT.index('GQ')
                AllSampleGQ=[ "." for i in range(0,len(SampleStrings))]
                for i in range(0,len(SampleStrings)):
                    if AllSampleGT[i]!="./.":
                        print AllSampleGT[i]
                        print SampleStrings[i]
                        print i
                        print linelist[0]+" "+linelist[1]
                        AllSampleGQ[i]=SampleStrings[i][GQind]
                        
            if "DP" in FORMAT:
            ## Define Allele Count
                DPind=FORMAT.index('DP')
                AllSampleDP=[ "." for i in range(0,len(SampleStrings))]
                for i in range(0,len(SampleStrings)):
                    if AllSampleGT[i]!="./.":
                        AllSampleDP[i]=SampleStrings[i][DPind]
            #headerlist=['Chromosome','Position','ID','REF','ALT','Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.ExAC','AlleleFrequency.1KG','AlleleFrequency.ESP','AlleleFrequency.VCF','SIFTprediction','PP2prediction','MAprediction','MTprediction','GERP++','CADDscore','MetaSVMprediction','SegmentalDuplication','PredictionSummary','VariantCallQuality','QUAL','MQ0','DP','FILTER']+[ i+" GT" for i in AllSamplesList]+[ i+" AD" for i in AllSamplesList]+['AlternateAlleles', 'INFO']
            OutputList=linelist[0:4]+[ALT,altnum+1,GeneName,VariantFunction,VariantClass,AAchange,ExACFreq,KGFreq,ESPFreq,VCFFreq,SIFTprediction,PP2prediction,MAprediction,MTprediction,GERPscore,CADDscore,MetaSVMprediction,SegDup,PathoLevel,FILTER,QUAL,MQ0number,DPnumber,VariantFilter]+AllSampleGT+AllSampleAD+AllSampleDP+AllSampleGQ+[AltAllsStr,INFOstring]
            OutputList= [ str(i) for i in OutputList ]
            OutputString="\t".join(OutputList)
            Output.write(OutputString+"\n")
            
Outlog.write("\t Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("Tabulation Complete: "+TimeNow+"\n")
print "Tabulation Complete "+TimeNow