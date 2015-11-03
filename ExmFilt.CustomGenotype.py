#!/usr/bin/env python
#$ -N CustFilt -cwd

################################################################################################################
##### USAGE                                                                                                  ###
################################################################################################################
#### The purpose of this script is to filter a VCF file for variants based on user-specified genotypes for each sample.
#### For example: To filter according to an autosomal dominant inheritance model, filter all mutations that are homozygous reference (0/0) in one parent and heterozygous (0/1) in the other parent and the proband.
#### It is not necessary to specify all of the samples in the vcf, if only a subset is supplied only that subset will be used for filtering and output to the tab delimited file; the output vcf will contain all samples. If no samples are specified then all variants will be returned
#### If a variant has multiple alternate alleles then the script will, with some limitations, iterate across possible combinations of these for heterozygous and homozygous alternate. Limitations:
####            The homozygous alternate allele will be the same for all individuals (i.e. if two individual are 2/2 and 3/3, this will not be kept)
####            The heterozygous genotype will always contain the homozyogus alternate genotype allele (i.e. if the homozygous atlernate is 3/3, a variant with the heterozygote 1/2 will not be kept but one with 1/3 or 2/3 will)
#### For example to run the autosomal dominant filter described above with the default quality filters:
####    ExmFilt.CustomFiltering.py -v MyData.vcf -o MyOutputFiles -e ProbandID,Parent1ID -r Parent2ID
################################################################################################################
################################################################################################################
################################################################################################################




################################################################################################################
##### GET PARAMETERS                                                                                         ###
################################################################################################################
import gzip
import os
from optparse import OptionParser
parser = OptionParser()

## Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="Input VCF file", metavar="VCFfile")
parser.add_option("-o", "--out", dest="OutputFileName",help="Base of names for output files", metavar="OutputFileName")

## Assign samples to heterozygous, homozygous or reference for genotype filtering
parser.add_option("-a", "--alt", dest="homozygousgenotypes",help="Samples that should be homozygous for alternate allele, e.g. 1/1, 2/2, separated by commas", metavar="homozygousgenotypes")
parser.add_option("-e", "--het", dest="heterozygousgenotypes",help="Samples that should be heterozygous, e.g. 0/1, separated by commas", metavar="heterozygousgenotypes")
parser.add_option("-r", "--ref", dest="referencegenotypes",help="Samples that should be homozygous for reference allele, i.e. 0/0, separated by commas", metavar="referencegenotypes")
parser.add_option("-x", "--xalt", dest="notalternate", help="Samples that should not be homozygous for alternate allele, i.e. maybe 0/0 or 0/1, separated by commas", metavar="notalternate")
parser.add_option("-n", "--xref", dest="notreference",help="Samples that should not be homozygous for reference allele, i.e. maybe  0/1 or 1/1, separated by commas", metavar="notreference")
parser.add_option("-u", "--unfl", dest="notfiltered",help="Samples that should not be used for filtereing but output the genoytpe data", metavar="notfiltered")

## Additional Parameters to Adjust filtering thresholds
parser.set_defaults(MQ0threshold="0.05")
parser.set_defaults(DPthreshold="3")
parser.set_defaults(ExACthreshold="0.01")
parser.set_defaults(OneKGthreshold="0.01")
parser.set_defaults(AllMAF="2")
parser.set_defaults(ESPthreshold="0.01")
parser.set_defaults(CHTthreshold="0.05")
parser.set_defaults(GQthreshold="30")
parser.set_defaults(HetAllCountthreshold="3")
parser.set_defaults(HomAllFracthreshold="0.7")


#### - unused options: b       l    q  t   w  z 

#arguments
parser.add_option("-d", "--dp", dest="DPthreshold",help="minimum depth of coverage", metavar="DPthreshold")
parser.add_option("-m", "--mq", dest="MQ0threshold",help="maximum for MQ value, given as a decimal fraction (e.g. 0.05 = 5% of reads for a variant can have MQ0) - fraction of reads with quality 0", metavar="MQ0")
parser.add_option("-f", "--maf", dest="AllMAF",help="set all Maximum MAF filters to this value", metavar="AllMAF")
parser.add_option("-y", "--exac", dest="ExACthreshold",help="maximum MAF in Exome Aggregation Consortium database", metavar="ExACthreshold")
parser.add_option("-k", "--okg", dest="OneKGthreshold",help="maximum MAF in One Thousand Genomes database", metavar="OneKGthreshold")
parser.add_option("-s", "--esp", dest="ESPthreshold",help="maximum MAF in ESP-GO database", metavar="ESPthreshold")
parser.add_option("-c", "--cht", dest="CHTthreshold",help="maximum frequency of allele in cohort of samples in the vcf", metavar="CHTthreshold")
parser.add_option("-g", "--gq", dest="GQthreshold",help="minimum for GQ value", metavar="GQthreshold")
parser.add_option("-i", "--aac", dest="HetAllCountthreshold",help="minimum alternate allele count in heterozygous calls", metavar="HetAllCountthreshold")
parser.add_option("-j", "--haf", dest="HomAllFracthreshold",help="Fraction of allele counts for homozygous allele", metavar="HomAllFracthreshold")

#Flags
parser.add_option("-P", "--nopathogenicity", action='store_true', dest="nopatho", help="Do not filter using pathogenicity predictions")
parser.add_option("-Q", "--noqual", action='store_true', dest="noqual", help="Do NOT filter by vcf FILTER field")

parser.add_option("-X", "--xlinked", action='store_true', dest="xlink", help="Only output X chromosome")
parser.add_option("-D", "--denovo", action='store_true', dest="denovo", help="Use higher stringency filters for de novo")
parser.add_option("-V", "--vcfout", action='store_true', dest="vcfout", help="Output the selected variants in vcf format")
parser.add_option("-S", "--keepsyn", action='store_true', dest="keepsyn", help="Keep synonymous variants as well")


parser.add_option("-Z", "--debug", action='store_true', dest="debug", help="Output a boolean decision file for debugging purposes")

(options, args) = parser.parse_args()

################################################################################################################
##### OPEN INPUT FILES AND OUTPUT FILES                                                                      ###
################################################################################################################

##Assign input and output files
FilTyp=str(options.VCFfile)
FilTyp=FilTyp.split(".")
FilTyp=str(FilTyp[-1])
if FilTyp == 'gz':
    VCF=gzip.open(options.VCFfile,'r')
elif FilTyp == 'vcf':
    VCF=open(options.VCFfile,'r')
else:
    print "Incorrect vcf file type"
    sys.exit(1)

BaseName=str(options.OutputFileName)

NoPatho=options.nopatho
NoQualityFilter=options.noqual
XLink=options.xlink
DeNovo=options.denovo
OutputVcf=options.vcfout
DeBug=options.debug
KeepSyn=options.keepsyn
TabOutputFilename=BaseName+'.tsv'
VcfOutputFilename=BaseName+'.vcf'
LogOutputFilename=BaseName+'.log'
PassOutputFilename=BaseName+'.boolean.log'
Output=open(TabOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')
if OutputVcf: Outvcf=open(VcfOutputFilename,'w')
if DeBug: OutPass=open(PassOutputFilename,'w')

################################################################################################################
##### ASSIGN FILTER VARIABLES                                                                                ###
################################################################################################################

## Variant Quality Filter

BadSnpFilters=['QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','LowQuality']
BadInDFilters=['LowQD_Indel','LowQuality']
MidSnpFilters=['FS_Mid_SNP','QD_Mid_SNP','VQSRTrancheSNP99.00to99.90','VQSRTrancheSNP99.90to100.00']
MidInDFilters=['FSBias_Indel','RPBias_Indel','VQSRTrancheINDEL99.00to99.90','VQSRTrancheINDEL99.90to100.00']

HetAllCountFilter=int(options.HetAllCountthreshold)
HomAllFrac=float(options.HomAllFracthreshold)
MQ0filter=float(options.MQ0threshold)
AllMAFfilters=float(options.AllMAF)
ExACfilter=float(options.ExACthreshold)
OneKGfilter=float(options.OneKGthreshold)
ESPfilter=float(options.ESPthreshold)
DPfilter=int(options.DPthreshold)
GQfilter=int(options.GQthreshold)
CHTfilter=float(options.CHTthreshold)
QualFilter=30
if AllMAFfilters != 2:
    ExACfilter=AllMAFfilters
    OneKGfilter=AllMAFfilters
    ESPfilter=AllMAFfilters

## Change filters if de novo required
if DeNovo:
    HetAllCountFilter=6
    HomAllFrac=0.98
    ExACfilter= ExACfilter / 10
    OneKGfilter= OneKGfilter / 10
    ESPfilter= ESPfilter / 10
    CHTfilter=0.01
    DPfilter=5

## Variant classes
MissenseClass=['nonsynonymousSNV','unknown']
NonFrameshiftClass=['nonframeshiftdeletion','nonframeshiftinsertion','nonframeshiftsubstitution']
FrameshiftClass=['frameshiftdeletion','frameshiftinsertion','frameshiftsubstitution']
InDelClass=NonFrameshiftClass+FrameshiftClass
NonsenseClass=['stopgain','stoploss']
SplicingClass=['splicing','exonic,splicing']

## Genotype Filters
AlternateSamples=options.homozygousgenotypes
HeterozygousSamples=options.heterozygousgenotypes
ReferenceSamples=options.referencegenotypes
NotAlternateSamples=options.notalternate
NotReferenceSamples=options.notreference
NotFilteredSamples=options.notfiltered
# initialise empty lists so that if no samples are specified we have an empty list for comparing to later
AlternateSampleList=[]
HeterozygousSampleList=[]
ReferenceSampleList=[]
NotAlternateSampleList=[]
NotReferenceSampleList=[]
NotFilteredSampleList=[]
if AlternateSamples is not None:
    AlternateSampleList=AlternateSamples.upper().split(',')
if HeterozygousSamples is not None:
    HeterozygousSampleList=HeterozygousSamples.upper().split(',')
if ReferenceSamples is not None:
    ReferenceSampleList=ReferenceSamples.upper().split(',')
if NotReferenceSamples is not None:
    NotReferenceSampleList=NotReferenceSamples.upper().split(',')
if NotAlternateSamples is not None:
    NotAlternateSampleList=NotAlternateSamples.upper().split(',')
if NotFilteredSamples is not None:
    NotFilteredSampleList=NotFilteredSamples.upper().split(',')

################################################################################################################
##### INITIALISE OUTPUT FILES                                                                                ###
################################################################################################################

##start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering log: "+TimeNow+"\n")
Outlog.write("VCF: "+str(options.VCFfile)+"\n")
Outlog.write("OutputName: "+str(options.OutputFileName)+"\n")
Outlog.write("\n")
Outlog.write("Genotype Filters: \n")
Outlog.write("\t Homozygous Alternate Samples: "+str(AlternateSamples)+"\n")
Outlog.write("\t Heterozygous Samples: "+str(HeterozygousSamples)+"\n")
Outlog.write("\t Homozygous Reference Samples: "+str(ReferenceSamples)+"\n")
Outlog.write("\t Not-Reference Samples: "+str(NotReferenceSamples)+"\n")
Outlog.write("\t Not-Alternate Samples: "+str(NotAlternateSamples)+"\n")
Outlog.write("\t Unfiltered Samples: "+str(NotFilteredSamples)+"\n")
Outlog.write("\n")

##Start output tables
AllSamplesList=AlternateSampleList+HeterozygousSampleList+ReferenceSampleList+NotReferenceSampleList+NotAlternateSampleList+NotFilteredSampleList
headerlist=['Chromosome','Position','ID','REF','ALT','Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.ExAC','AlleleFrequency.1KG','AlleleFrequency.ESP','MetaSVM','SIFTprediction','PP2prediction','MAprediction','MTprediction','GERP++','CADDscore','SegmentalDuplication','Cosmic','PredictionSummary','VariantCallQuality']+[ i+" GT" for i in AllSamplesList]+[ i+" AD" for i in AllSamplesList]+['AlternateAlleles', 'MetaSVMScore','FILTER', 'INFO']+AllSamplesList+['ESP.aa','ESP.ea','1KG.eur','1KG.amr','1KG.eas','1KG.afr','1KG.sas','ExAC.afr','ExAC.amr','ExAC.eas','ExAC.fin','ExAC.nfe','ExAC.oth','ExAC.sas']
Output.write("\t".join(headerlist)+"\n")
if DeBug:
    debugheaderlist=['Chromosome','Position','REF','ALT','AAlt','PassExAC','PassKG','PassESP','PassVCF','PassQUAL','PassMQ','PassFunction','PassClass','PassGT','PassDP','PassGQ','PassPatho','PassFilter','PassAFC']
    OutPass.write("\t".join(debugheaderlist)+"\n")
    

################################################################################################################
##### INITIALISE COUNTS OF VARIANTS AND START FILTERING                                                      ###
################################################################################################################

## Read VCF file
NameToColumn={}
ColumnToName={}
OrigCount=0
FiltCount=[0,0,0]
SilentCount=[0,0,0]
MissenseCount=[0,0,0]
NonFrameshiftCount=[0,0,0]
FrameshiftCount=[0,0,0]
NonsenseCount=[0,0,0]
SplicingCount=[0,0,0]
UnknownCount=[0,0,0]

for line in VCF:
    if '#' in line and OutputVcf: Outvcf.write(line)
    if  line.startswith('#CHROM'):
    ################################################################################################################
    ##### MAP COLUMN NAME TO NUMBER AND THEN LIST COLUMN NUMBERS CORRESPONDING TO EACH SAMPLE                    ###
    ################################################################################################################
        ColumnAndNumber=enumerate(line.strip().upper().split("\t"))
        for i in ColumnAndNumber:
            NameToColumn[i[1]]=i[0]
            ColumnToName[i[0]]=i[1]
        AlternateColumnNumber=[ NameToColumn[i] for i in AlternateSampleList ]
        HeterozygousColumnNumber=[ NameToColumn[i] for i in HeterozygousSampleList ]
        ReferenceColumnNumber=[ NameToColumn[i] for i in ReferenceSampleList ]
        NotReferenceColumnNumber=[ NameToColumn[i] for i in NotReferenceSampleList ]
        NotAlternateColumnNumber=[ NameToColumn[i] for i in NotAlternateSampleList ]
        NotFilteredColumnNumber=[ NameToColumn[i] for i in NotFilteredSampleList ]
        #Adjust Cohort Allele Frequencies to compensate for low sample numbers:
        linelist=line.split("\t")
        SampleNum=len(linelist)-9
        AltCnt= len(AlternateSampleList) * 2
        HetCnt= len(HeterozygousSampleList)
        FamAllCnt= HetCnt + ( 2 * AltCnt )
        TotAllNum= SampleNum * 2
        ExtrAllCnt= (TotAllNum - FamAllCnt) * CHTfilter
        TotAllCnt= FamAllCnt + ExtrAllCnt + 0.5 # plus 0.5 to compensate for rounding errors and ensure always at least 1 variant is allowed
        minMAF=round(TotAllCnt/TotAllNum,2)
        CHTfilter=max(CHTfilter, minMAF)
    if '#' not in line:
        ## Skip all chromosomes except X is X-linked is specified
        linelist=line.split("\t", 1)
        Chrom=linelist[0]
        PassX=True
        if Chrom!="X" and XLink:
            PassX=False
    if '#' not in line and PassX:
        ################################################################################################################
        ##### PARSE LINE AND GET VARIANT INFO                                                                        ###
        ################################################################################################################
        
        OrigCount=OrigCount+1
        linelist=line.split("\t")
        
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
        
        #Variant details
        GeneName=INFOdict.get('GeneName','.')
        VariantFunction=INFOdict.get('VarFunc','none')
        VariantClassList=INFOdict.get('VarClass','none').split(',')
        AAchangeList=INFOdict.get('AAChange','.').split(',')
        
        #Calling quality
        DPnumber=float(INFOdict.get('DP','1'))
        MQ0number=float(INFOdict.get('MQ0','0'))
        SegDup=INFOdict.get('SegDup','none')
        
        #Allele Frequencies
        KGFreqList=str(INFOdict.get('1KGfreq','.')).split(',')
        ESPFreqList=str(INFOdict.get('ESPfreq','.')).split(',')
        ExACFreqList=str(INFOdict.get('ExACfreq','.')).split(",")
        VCFFreqList=str(INFOdict.get('AF',0)).split(",")
        
        #population specific
        ESPaaFreqList=str(INFOdict.get('ESP.aa.freq','.')).split(',')
        ESPeaFreqList=str(INFOdict.get('ESP.ea.freq','.')).split(',')
        KGeurFreqList=str(INFOdict.get('1KG.eur.freq','.')).split(',')
        KGamrFreqList=str(INFOdict.get('1KG.amr.freq','.')).split(',')
        KGeasFreqList=str(INFOdict.get('1KG.eas.freq','.')).split(',')
        KGafrFreqList=str(INFOdict.get('1KG.afr.freq','.')).split(',')
        KGsasFreqList=str(INFOdict.get('1KG.sas.freq','.')).split(',')
        ExACafrFreqList=str(INFOdict.get('ExAC.afr.freq','.')).split(',')
        ExACamrFreqList=str(INFOdict.get('ExAC.amr.freq','.')).split(',')
        ExACeasFreqList=str(INFOdict.get('ExAC.eas.freq','.')).split(',')
        ExACfinFreqList=str(INFOdict.get('ExAC.fin.freq','.')).split(',')
        ExACnfeFreqList=str(INFOdict.get('ExAC.nfe.freq','.')).split(',')
        ExACothFreqList=str(INFOdict.get('ExAC.oth.freq','.')).split(',')
        ExACsasFreqList=str(INFOdict.get('ExAC.sas.freq','.')).split(',')
        
        #Variant effect predictors
        Cosmic=INFOdict.get('COSMIC','none')
        SIFTpredictionList=INFOdict.get('SIFTprd','.').split(',')
        PP2predictionList=INFOdict.get('PP2.hvar.prd','.')
        PP2predictionList=PP2predictionList.split(',')
        MApredictionList=INFOdict.get('MutAprd','.').split(',')
        MTpredictionList=INFOdict.get('MutTprd','.').split(',')
        GERPscoreList=str(INFOdict.get('GERP','.')).split(',')
        CADDscoreList=str(INFOdict.get('CADDphred','.')).split(',')
        SVMscoreList=INFOdict.get('MetaSVMscr','.').split(',')
        SVMpredictionList=INFOdict.get('MetaSVMprd','.').split(',')
        
        ## Get Sample Level Variant information
        AlternateQualityString=[ linelist[i].strip() for i in AlternateColumnNumber ]
        HeterozygousQualityString=[ linelist[i].strip() for i in HeterozygousColumnNumber ]
        ReferenceQualityString=[ linelist[i].strip() for i in ReferenceColumnNumber ]
        NotReferenceQualityString=[ linelist[i].strip() for i in NotReferenceColumnNumber ]
        NotAlternateQualityString=[ linelist[i].strip() for i in NotAlternateColumnNumber ]
        NotFilteredQualityString=[ linelist[i].strip() for i in NotFilteredColumnNumber ]
        
        AlternateQualityList=[ i.split(':') for i in AlternateQualityString ]
        HeterozygousQualityList=[ i.split(':') for i in HeterozygousQualityString ]
        ReferenceQualityList=[ i.split(':') for i in ReferenceQualityString ]
        NotReferenceQualityList=[ i.split(':') for i in NotReferenceQualityString ]
        NotAlternateQualityList=[ i.split(':') for i in NotAlternateQualityString ]
        NotFilteredQualityList=[ i.split(':') for i in NotFilteredQualityString ]
        
        ## Define Genotypes
        AlternateGT=[ AlternateQualityList[i][0] for i in range(0,len(AlternateQualityString)) ]
        HeterozygousGT=[ HeterozygousQualityList[i][0] for i in range(0,len(HeterozygousQualityString)) ]
        ReferenceGT=[ ReferenceQualityList[i][0] for i in range(0,len(ReferenceQualityString)) ]
        NotReferenceGT=[ NotReferenceQualityList[i][0] for i in range(0,len(NotReferenceQualityString)) ]
        NotAlternateGT=[ NotAlternateQualityList[i][0] for i in range(0,len(NotAlternateQualityString)) ]
        NotFilteredGT=[ NotFilteredQualityList[i][0] for i in range(0,len(NotFilteredQualityString)) ]
        AllSampleGT=AlternateGT+HeterozygousGT+ReferenceGT+NotReferenceGT+NotAlternateGT+NotFilteredGT
        
        ################################################################################################################
        ##### INITIAL LOCUS LEVEL FILTERS                                                                            ###
        ################################################################################################################

        ##Check Variant Quality
        
        PassQUAL=False
        if QUAL >= QualFilter:
            PassQUAL=True
        
        ## Check if MQ0 passes threshold
        PassMQ=False
        MQ0Fraction=float(MQ0number)/float(DPnumber)
        if MQ0Fraction <= MQ0filter:
            PassMQ=True
         
        ## Check if Variant class passes
        PassFunction=False
        if VariantFunction=='exonic' or VariantFunction in SplicingClass or VariantFunction=='none':
            PassFunction=True
        
        ##Check FILTER field passes - separately for SNPs and Indels, will check which class is relevant later
        GoodSNP=False
        GoodInDel=False
        if all( str(i) not in VariantFilter for i in BadInDFilters):
            GoodInDel=True
        if all( str(i) not in VariantFilter for i in BadSnpFilters):
            GoodSNP=True
        
        if NoQualityFilter:
            PassMQ=True
            PassQUAL=True
            GoodInDel=True
            GoodSNP=True
        
        ################################################################################################################
        ##### SAMPLE AND VARIANT LEVEL FILTERS                                                                       ###
        ################################################################################################################
        
        ## Check if all genotypes are present
        PassGeno=False
        if './.' not in AlternateGT and './.' not in HeterozygousGT and './.' not in ReferenceGT and './.' not in NotReferenceGT and './.' not in NotAlternateGT:
            PassGeno=True
        
        if PassGeno and PassFunction and PassMQ and PassQUAL:
            ################################################################################################################
            ##### SAMPLES - DP and GQ                                                                                    ###
            ################################################################################################################

            PassDP=True
            if "DP" in FORMAT:
                PassDP=False
                ## Define DP
                DPind=FORMAT.index('DP')
                AlternateDP=[ str(AlternateQualityList[i][DPind]) for i in range(0,len(AlternateQualityString)) ]
                for i in range(0,len(AlternateDP)):
                    if AlternateDP[i]==".":
                        AlternateDP[i]="0"
                    AlternateDP[i]=float(AlternateDP[i])
                HeterozygousDP=[ str(HeterozygousQualityList[i][DPind]) for i in range(0,len(HeterozygousQualityString)) ]
                for i in range(0,len(HeterozygousDP)):
                    if HeterozygousDP[i]==".":
                        HeterozygousDP[i]="0"
                    HeterozygousDP[i]=float(HeterozygousDP[i])
                ReferenceDP=[ str(ReferenceQualityList[i][DPind]) for i in range(0,len(ReferenceQualityString)) ]
                for i in range(0,len(ReferenceDP)):
                    if ReferenceDP[i]==".":
                        ReferenceDP[i]="0"
                    ReferenceDP[i]=float(ReferenceDP[i])
                NotReferenceDP=[ str(NotReferenceQualityList[i][DPind]) for i in range(0,len(NotReferenceQualityString)) ]
                for i in range(0,len(NotReferenceDP)):
                    if NotReferenceDP[i]==".":
                        NotReferenceDP[i]="0"
                    NotReferenceDP[i]=float(NotReferenceDP[i])
                NotAlternateDP=[ str(NotAlternateQualityList[i][DPind]) for i in range(0,len(NotAlternateQualityString)) ]
                for i in range(0,len(NotAlternateDP)):
                    if NotAlternateDP[i]==".":
                        NotAlternateDP[i]="0"
                    NotAlternateDP[i]=float(NotAlternateDP[i])
                ## Check to see if Depth passes
                if all(i >= DPfilter for i in AlternateDP) and all(i >= DPfilter for i in HeterozygousDP) and all(i >= DPfilter for i in ReferenceDP):
                    PassDP=True
            
            PassGQ=True
            if "GQ" in FORMAT:
                PassGQ=False
                ## Define GQ
                GQind=FORMAT.index('GQ')
                AlternateGQ=[ float(AlternateQualityList[i][GQind]) for i in range(0,len(AlternateQualityString)) ]
                HeterozygousGQ=[ float(HeterozygousQualityList[i][GQind]) for i in range(0,len(HeterozygousQualityString)) ]
                ReferenceGQ=[ float(ReferenceQualityList[i][GQind]) for i in range(0,len(ReferenceQualityString)) ]
                NotReferenceGQ=[ float(NotReferenceQualityList[i][GQind]) for i in range(0,len(NotReferenceQualityString)) ]
                NotAlternateGQ=[ float(NotAlternateQualityList[i][GQind]) for i in range(0,len(NotAlternateQualityString)) ]
                ## Check to see if genotype quality passes
                if all(i >= GQfilter for i in AlternateGQ) and all(i >= GQfilter for i in HeterozygousGQ) and all(i >= GQfilter for i in ReferenceGQ) and all(i >= GQfilter for i in NotReferenceGQ) and all(i >= GQfilter for i in NotAlternateGQ):
                    PassGQ=True
            ################################################################################################################
            ##### CHECK TO SEE IF GENOTYPES PASS, ITERATE ACROSS MULTIPLE ALT ALLELES IF NECESSARY                      ###
            ################################################################################################################
            
            AltRng=range(0, AltNum)
            for altnum in AltRng:
                
                ################################################################################################################
                ##### CHECK GENOTYPES                                                                                         ###
                ################################################################################################################
                
                ##check GT
                PassGT=False
                AltAll=str(altnum+1)
                if all(i.count(AltAll)==2 for i in AlternateGT) and all(i.count(AltAll)==1 for i in HeterozygousGT) and all(i.count(AltAll)==0 for i in ReferenceGT) and all(i.count(AltAll)>0 for i in NotReferenceGT) and all(i.count(AltAll)<2 for i in NotAlternateGT):
                    PassGT=True
                
            
                ################################################################################################################
                ##### PER SAMPLE ALTERNATE ALLELE COUNTS/FRACTIONS                                                           ###
                ################################################################################################################
                
                PassAFC=True
                AllSampleAD=[ "NA" for i in AllSampleGT] 
                if "AD" in FORMAT:
                    PassAFC=False
                    ## Define Allele Count
                    ADind=FORMAT.index('AD')
                    
                    AlternateAC=[ AlternateQualityList[i][ADind] for i in range(0,len(AlternateQualityString)) ]
                    AllSampleAD=AlternateAC
                    AlternateAC=[ i.split(',') for i in AlternateAC ]
                    AlternateTotal=[ sum(int(j) for j in i)  for i in AlternateAC ]
                    AlternateAC=[ int(i[altnum+1]) for i in AlternateAC ]
                    AlternateAAF=[ float(i)/float(max(j,1)) for i,j in zip(AlternateAC,AlternateTotal) ] ##if Total is 0 then python throws and error, so set to 0, this will still leave the AAF as 0 as if the total is 0 then the AC must also be 0
                    
                    HeterozygousAC=[ HeterozygousQualityList[i][ADind] for i in range(0,len(HeterozygousQualityString)) ]
                    AllSampleAD=AllSampleAD+HeterozygousAC
                    HeterozygousAC=[ i.split(',') for i in HeterozygousAC ]
                    HeterozygousAC=[ int(i[altnum+1]) for i in HeterozygousAC ]
                    
                    ReferenceAC=[ ReferenceQualityList[i][ADind] for i in range(0,len(ReferenceQualityString)) ]
                    AllSampleAD=AllSampleAD+ReferenceAC
                    ReferenceAC=[ i.split(',') for i in ReferenceAC ]
                    ReferenceTotal=[ sum(int(j) for j in i)  for i in ReferenceAC ]
                    ReferenceAC=[ int(i[0]) for i in ReferenceAC ]
                    ReferenceAAF=[ float(i)/float(max(j,1)) for i,j in zip(ReferenceAC,ReferenceTotal) ]
                    
                    NotReferenceAC=[ NotReferenceQualityList[i][ADind] for i in range(0,len(NotReferenceQualityString)) ]
                    AllSampleAD=AllSampleAD+NotReferenceAC
                    NotReferenceAC=[ i.split(',') for i in NotReferenceAC ]
                    NotReferenceTotal=[ sum(int(j) for j in i)  for i in NotReferenceAC ]
                    NotReferenceAC=[ int(i[altnum+1]) for i in NotReferenceAC ]
                    NotReferenceAAF=[ float(i)/float(max(j,1)) for i,j in zip(NotReferenceAC,NotReferenceTotal) ]
                    
                    NotAlternateAC=[ NotAlternateQualityList[i][ADind] for i in range(0,len(NotAlternateQualityString)) ]
                    AllSampleAD=AllSampleAD+NotAlternateAC
                    NotAlternateAC=[ i.split(',') for i in NotAlternateAC ]
                    NotAlternateTotal=[ sum(int(j) for j in i)  for i in NotAlternateAC ]
                    NotAlternateAAC=[ int(i[altnum+1]) for i in NotAlternateAC ]
                    NotAlternateRAC=[ int(i[0]) for i in NotAlternateAC ]
                    NotAlternateRAF=[ float(i)/float(max(j,1)) for i,j in zip(NotAlternateRAC,NotAlternateTotal) ]
                    
                    NotFilteredAC=list()
                    for i in range(0,len(NotFilteredQualityString)):
                        if './.' not in NotFilteredQualityList[i]:
                            NotFilteredAC.append(NotFilteredQualityList[i][ADind])
                        else:
                            NotFilteredAC.append(0)
                    AllSampleAD=AllSampleAD+NotFilteredAC
                    
                    ##Check alternate allele counts and fractions
                    if all( i >= HetAllCountFilter for i in HeterozygousAC) and all( i >= HomAllFrac for i in AlternateAAF) and all( i >= HomAllFrac for i in ReferenceAAF) and all( ( i.count(AltAll)==1 and j >= HetAllCountFilter ) or ( i.count(AltAll)==2 and k >= HomAllFrac ) for i,j,k in zip(NotReferenceGT,NotReferenceAC,NotReferenceAAF) ) and all( ( i.count(AltAll)==1 and j >= HetAllCountFilter ) or ( i.count(AltAll)==0 and k >= HomAllFrac ) for i,j,k in zip(NotAlternateGT,NotAlternateAC,NotAlternateRAF) ):
                        PassAFC=True
                    
                ################################################################################################################
                ##### VARIANT SPECIFIC CLASS CHECK, FREQUENCY AND VARIANT CALLING QUALITY                                    ###
                ################################################################################################################
                
                ##check variant class
                PassClass=False
                cltnum=min(len(VariantClassList)-1, altnum)
                VariantClass=str(VariantClassList[cltnum])
                if VariantClass != 'synonymousSNV':
                    PassClass=True
                
                ## a function to get the requisite frequency value, if the annotation was missing in the vcf then we have only a single "." so need to adjust the index number, also for "." return 0
                def getInfoVal ( InfoFieldList ):
                    cltnum=min(len(InfoFieldList)-1, altnum)
                    InfoField=str(InfoFieldList[cltnum])
                    if InfoField == ".":
                        InfoFieldTest=float(0)
                    else:
                        InfoFieldTest=float(InfoField)
                    return [InfoField,InfoFieldTest]
                
                
                
                ## Check if KG passes threshold
                PassKG=False
                KGFreq=getInfoVal(KGFreqList)[0]
                KGFreqtest=getInfoVal(KGFreqList)[1]
                if KGFreqtest <= OneKGfilter:
                    PassKG=True
                
                ## Check if ESP passes threshold
                PassESP=False
                ESPFreq=getInfoVal(ESPFreqList)[0]
                ESPFreqtest=getInfoVal(ESPFreqList)[1]
                if ESPFreqtest <= ESPfilter:
                    PassESP=True
                    
                ## Check if ExAC passes threshold
                PassExAC=False
                ExACFreq=getInfoVal(ExACFreqList)[0]
                ExACFreqtest=getInfoVal(ExACFreqList)[1]
                if ExACFreqtest <= ExACfilter:
                    PassExAC=True
                
                ## Check if VCF passes threshold
                PassVCF=False
                VCFFreqtest=getInfoVal(VCFFreqList)[1]
                if VCFFreqtest <= CHTfilter:
                    PassVCF=True
                
                #get the population specific frequencies
                ESPaaFreq=getInfoVal(ESPaaFreqList)[0]
                ESPeaFreq=getInfoVal(ESPeaFreqList)[0]
                KGeurFreq=getInfoVal(KGeurFreqList)[0]
                KGamrFreq=getInfoVal(KGamrFreqList)[0]
                KGeasFreq=getInfoVal(KGeasFreqList)[0]
                KGafrFreq=getInfoVal(KGafrFreqList)[0]
                KGsasFreq=getInfoVal(KGsasFreqList)[0]
                ExACafrFreq=getInfoVal(ExACafrFreqList)[0]
                ExACamrFreq=getInfoVal(ExACamrFreqList)[0]
                ExACeasFreq=getInfoVal(ExACeasFreqList)[0]
                ExACfinFreq=getInfoVal(ExACfinFreqList)[0]
                ExACnfeFreq=getInfoVal(ExACnfeFreqList)[0]
                ExACothFreq=getInfoVal(ExACothFreqList)[0]
                ExACsasFreq=getInfoVal(ExACsasFreqList)[0]

                
                ## Check and annotate pathogenicity (discard missense predicted to be benign - need to adjust this in splicining regions)
                
                cltnum=min(len(SIFTpredictionList)-1, altnum)
                SIFTprediction=SIFTpredictionList[cltnum]
                cltnum=min(len(PP2predictionList)-1, altnum)
                PP2prediction=PP2predictionList[cltnum]
                cltnum=min(len(CADDscoreList)-1, altnum)
                CADDscore=str(CADDscoreList[cltnum])
                cltnum=min(len(SVMpredictionList)-1, altnum)
                SVMprediction=SVMpredictionList[cltnum]
                if CADDscore == ".":
                    CADDscoretest=float(0)
                else:
                    CADDscoretest=float(CADDscore)
                
                PassPatho=False
                PathoLevel="Low"
                if VariantClass in MissenseClass and SIFTprediction=="." and PP2prediction=="." and CADDscore=="." and SVMprediction==".":
                    PathoLevel="Med"
                    PassPatho=True
                if VariantClass in MissenseClass and (SIFTprediction=="D" or PP2prediction=="D" or PP2prediction=="P" or CADDscoretest>=15):
                    PathoLevel="Med"
                    PassPatho=True
                if VariantClass in MissenseClass and (SIFTprediction=="D" and PP2prediction=="D"):
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass in MissenseClass and (SVMprediction=="D"):
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass in MissenseClass and (CADDscoretest>=25):
                    PathoLevel="High"
                    PassPatho=True
                if VariantFunction in SplicingClass:
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass in InDelClass or VariantClass in NonsenseClass:
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass == "synonymousSNV":
                    PathoLevel="Silent"
                
                ##check variant specific FILTER
                PassFilter=False
                if VariantClass in InDelClass and GoodInDel:
                    PassFilter=True
                if VariantClass not in InDelClass and GoodSNP:
                    PassFilter=True
                    
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
                
                ##Filter Overrides
                if NoQualityFilter:
                    PassDP=True
                    PassGQ=True
                    PassFilter=True
                    PassAFC=True
                if NoPatho:
                    PassPatho=True
                if KeepSyn and VariantClass == "synonymousSNV":
                    PassClass=True
                    PassPatho=True
                    
                ################################################################################################################
                ##### OUTPUTS                                                                                                ###
                ################################################################################################################
                
                ## Output Debug file if requested
                if DeBug:
                    OutPass.write("\t"+linelist[0]+" "+linelist[1]+" "+linelist[3]+" "+linelist[4]+" "+str(altnum)+" "+str(PassExAC)+" "+str(PassKG)+" "+str(PassESP)+" "+str(PassVCF)+" "+str(PassQUAL)+" "+str(PassMQ)+" "+str(PassFunction)+" "+str(PassClass)+" "+str(PassGT)+" "+str(PassDP)+" "+str(PassGQ)+" "+str(PassPatho)+" "+str(PassFilter)+" "+str(PassAFC)+"\n")
                ## If all pass then output line
                if PassKG and PassESP and PassExAC and PassVCF and PassClass and PassGT and PassDP and PassGQ and PassFilter and PassAFC and PassPatho:
                    cltnum=min(len(AAchangeList)-1, altnum)
                    AAchange=AAchangeList[cltnum]
                    cltnum=min(len(MApredictionList)-1, altnum)
                    MAprediction=MApredictionList[cltnum]
                    cltnum=min(len(MTpredictionList)-1, altnum)
                    MTprediction=MTpredictionList[cltnum]
                    cltnum=min(len(SVMscoreList)-1, altnum)
                    SVMscore=SVMscoreList[cltnum]                    
                    ##output
                    ALT=str(AltAlls[altnum])
                    cltnum=min(len(GERPscoreList)-1, altnum)
                    GERPscore=str(GERPscoreList[cltnum])
                    OutputList=linelist[0:4]+[ALT,GeneName,VariantFunction,VariantClass,AAchange,ExACFreq,KGFreq,ESPFreq,SVMprediction,SIFTprediction,PP2prediction,MAprediction,MTprediction,GERPscore,CADDscore,SegDup,Cosmic,PathoLevel,FILTER]+AllSampleGT+AllSampleAD+[AltAllsStr,SVMscore,VariantFilter,INFOstring]+AlternateQualityString+HeterozygousQualityString+ReferenceQualityString+NotReferenceQualityString+NotAlternateQualityString+NotFilteredQualityString+[ESPaaFreq,ESPeaFreq,KGeurFreq,KGamrFreq,KGeasFreq,KGafrFreq,KGsasFreq,ExACafrFreq,ExACamrFreq,ExACeasFreq,ExACfinFreq,ExACnfeFreq,ExACothFreq,ExACsasFreq]
                    OutputList= [ str(i) for i in OutputList ]
                    OutputString="\t".join(OutputList)
                    Output.write(OutputString+"\n")
                    if OutputVcf: Outvcf.write(line)
                    
                    ################################################################################################################
                    ##### ADD TO VARIANT COUNTERS                                                                                  ###
                    ################################################################################################################
                    
                    cntnum=0
                    if PathoLevel=="Medium":
                        cntnum=1
                    if PathoLevel=="High":
                        cntnum=2
                    FiltCount[cntnum] = FiltCount[cntnum] + 1
                    if VariantClass == "synonymousSNV":
                        SilentCount[cntnum] = SilentCount[cntnum] + 1
                    if VariantClass in MissenseClass:
                        MissenseCount[cntnum] = MissenseCount[cntnum] + 1
                    if VariantClass in NonFrameshiftClass:
                        NonFrameshiftCount[cntnum] = NonFrameshiftCount[cntnum] + 1
                    if VariantClass in FrameshiftClass:
                        FrameshiftCount[cntnum] = FrameshiftCount[cntnum] + 1
                    if VariantClass in NonsenseClass:
                        NonsenseCount[cntnum] = NonsenseCount[cntnum] + 1
                    if VariantFunction in SplicingClass:
                        SplicingCount[cntnum] = SplicingCount[cntnum] + 1
                    if VariantClass=="Unknown":
                        UnknownCount[cntnum] = UnknownCount[cntnum] + 1
################################################################################################################
##### OUTPUT Filter details to Log                                                                           ###
################################################################################################################
if not NoQualityFilter:
    Outlog.write("Individual Sample Filters: \n")
    Outlog.write("\t Minimum Depth of Coverage: "+str(DPfilter)+"\n")
    Outlog.write("\t Minimum Genotyping Quality(GQ): "+str(GQfilter)+"\n")
    Outlog.write("\t Minimum Heterozygous allele count: "+str(HetAllCountFilter)+"\n")
    Outlog.write("\t Minimum Homozygous allele fraction: "+str(HomAllFrac)+"\n")
Outlog.write("Variant Filters: \n")
if not NoQualityFilter:
    Outlog.write("\t MQ0/DP maximum: "+str(MQ0filter)+"\n")
if NoQualityFilter:
    Outlog.write("\t No Variant Quality Filters\n")
Outlog.write("\t 1000 genomes alternate allele frequency maximum: "+str(OneKGfilter)+"\n")
Outlog.write("\t GO ESP alternate allele frequency maximum: "+str(ESPfilter)+"\n")
Outlog.write("\t ExAC alternate allele frequency maximum: "+str(ExACfilter)+"\n")
Outlog.write("\t Within VCF allele frequency maximum: "+str(CHTfilter)+"\n")
if NoPatho:
    Outlog.write("\t No Predicted Pathogenecity Filters")
if not NoPatho:
    Outlog.write("\t Predicted Pathogenecity Filters: at least SIFT=D or PP2=D/P or CADD >= 15; all InDels and Nonsense \n")
Outlog.write("\n")
Outlog.write("\n")                        
################################################################################################################
##### OUTPUT VARIANT COUNTS TO LOG                                                                           ###
################################################################################################################
Outlog.write("---------------------------------------------------------------------\n")
Outlog.write("Results: \n")
Outlog.write("\t Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("\t Number of variants selected (Low/Med/High/Total): "+str(FiltCount[0])+"/"+str(FiltCount[1])+"/"+str(FiltCount[2])+"/"+str(sum(FiltCount))+"\n")
Outlog.write("\t Number of Missense variants selected (Med/High/Total): "+str(MissenseCount[0])+"/"+str(MissenseCount[1])+"/"+str(MissenseCount[2])+"/"+str(sum(MissenseCount))+"\n")
Outlog.write("\t Number of Nonsense variants selected (Med/High/Total): "+str(NonsenseCount[0])+"/"+str(NonsenseCount[1])+"/"+str(NonsenseCount[2])+"/"+str(sum(NonsenseCount))+"\n")
Outlog.write("\t Number of NonFrameshift InDels selected (Med/High/Total): "+str(NonFrameshiftCount[0])+"/"+str(NonFrameshiftCount[1])+"/"+str(NonFrameshiftCount[2])+"/"+str(sum(NonFrameshiftCount))+"\n")
Outlog.write("\t Number of Frameshift InDels selected (Med/High/Total): "+str(FrameshiftCount[0])+"/"+str(FrameshiftCount[1])+"/"+str(FrameshiftCount[2])+"/"+str(sum(FrameshiftCount))+"\n")
Outlog.write("\t Number of Splice site variants selected (Med/High/Total): "+str(SplicingCount[0])+"/"+str(SplicingCount[1])+"/"+str(SplicingCount[2])+"/"+str(sum(SplicingCount))+"\n")
Outlog.write("\t Number of Silent function variants selected (Med/High/Total): "+str(SilentCount[0])+"/"+str(SilentCount[1])+"/"+str(SilentCount[2])+"/"+str(sum(SilentCount))+"\n")
Outlog.write("\t Number of Unknown function variants selected (Med/High/Total): "+str(UnknownCount[0])+"/"+str(UnknownCount[1])+"/"+str(UnknownCount[2])+"/"+str(sum(UnknownCount))+"\n")
Outlog.write("\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering Finished: "+TimeNow+"\n")
print "Filtering finished "+TimeNow
