#!/usr/bin/env python
#$ -N CustFilt -cwd

################################################################################################################
##### USAGE                                                                                                  ###
################################################################################################################
#### The purpose of this script is to filter a VCF file for all variants in a specific gene found in a group of samples
#### For example: To filter according to an autosomal dominant inheritance model, filter all mutations that are homozygous reference (0/0) in one parent and heterozygous (0/1) in the other parent and the proband.
####    -v/--vcf    <required>    The script requires an input vcf
####    -o/--out    <required>    The use should specify a base name for the output files.The script outputs the filtered results as both a tab delimited file and a vcf file. The script also produces a log file and an auxilary file that contains the results of every test of every variant as "TRUE" or "FALSE" (primarily for debugging purposes).
####    -s/--sam    Sample ids as in the vcf separated by a comma. Leave blank for all variants and no sample data
####    -g/--gene   <required>    Gene Symbol to search for in the INFO field
####
################################################################################################################
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

## Assign samples to heterozygous, homozygous or reference for genotype filtering
parser.set_defaults(sampleIDs="ALL")
parser.add_option("-s", "--sam", dest="sampleIDs",help="Samples IDs, must match vcf CHROM line", metavar="sampleIDs")
parser.add_option("-g", "--gene", dest="geneName",help="Gene Symbol to search for in the INFO field", metavar="geneName")

(options, args) = parser.parse_args()

################################################################################################################
##### OPEN INPUT FILES AND INITIALISE OUTPUT FILES                                                           ###
################################################################################################################

##Assign input and output files
VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
TabOutputFilename=BaseName+'.tsv'
VcfOutputFilename=BaseName+'.vcf'
LogOutputFilename=BaseName+'.log'
Output=open(TabOutputFilename,'w')
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')

################################################################################################################
##### ASSIGN FILTER VARIABLES                                                                                ###
################################################################################################################

## Gene/Text
SearchGene=str(options.geneName)
# initialise empty lists so that if no samples are specified we have an empty list for comparing to later
Samples=options.sampleIDs
Sample=[]
if Samples is not None:
    Sample=Samples.split(',') #.upper().split(',')

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
Outlog.write("Samples: "+str(Samples)+"\n")
Outlog.write("Gene: "+str(SearchGene)+"\n")

#adjust search term to get exact matchs only
SearchGene='='+SearchGene+";"

#initialise arrays
NameToColumn={}
ColumnToName={}
BadGT=['./.', '0/0']

for line in VCF:
    if '#' in line:
        Outvcf.write(line)
    if  line.startswith('#CHROM'):
    ################################################################################################################
    ##### MAP COLUMN NAME TO NUMBER AND THEN LIST COLUMN NUMBERS CORRESPONDING TO EACH SAMPLE                    ###
    ################################################################################################################
        linelist=line.strip().split("\t")
        ColumnAndNumber=enumerate(linelist)
        print Sample
        if Sample[0]=="ALL":
            Sample=linelist[9:]
            print Sample
        for i in ColumnAndNumber:
            NameToColumn[i[1]]=i[0]
            ColumnToName[i[0]]=i[1]
        ColumnNumber=[ NameToColumn[i] for i in Sample ]
        ##Start output tables
        headerlist=['Chromosome','Position','ID','REF','ALT','Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.ExAC','AlleleFrequency.1KG','AlleleFrequency.ESP','SIFTprediction','PP2prediction','MAprediction','MTprediction','GERP++','CADDscore','SegmentalDuplication']+Sample+['FILTER', 'INFO']+Sample
        Output.write("\t".join(headerlist)+"\n")
    if '#' not in line:
        ################################################################################################################
        ##### PARSE LINE AND GET VARIANT INFO                                                                        ###
        ################################################################################################################
        
        linelist=line.split("\t")
        
        ##Get Alternate Alleles
        AltAlls=linelist[4]
        AltAlls=AltAlls.split(",")
        AltNum=len(AltAlls)
        
        QUAL=float(linelist[5])
        FILTER=linelist[6]
        
        ## Get variant data from INFO Field
        INFOstring=linelist[7]
        INFOcolumn=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumn:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        DPnumber=float(INFOdict.get('DP','.'))
        MQ0number=float(INFOdict.get('MQ0','.'))
        GeneName=INFOdict.get('GeneName','.')
        VariantFunction=INFOdict.get('VarFunc','none')
        SegDup=INFOdict.get('SegDup','none')
        
        KGFreq=str(INFOdict.get('1KGfreq','.'))
        ESPFreq=str(INFOdict.get('ESPfreq','.'))
        ExACFreq=str(INFOdict.get('ExACfreq','.'))
        VCFFreq=str(INFOdict.get('AF',0))
        
        VariantClass=INFOdict.get('VarClass','none')
        AAchange=INFOdict.get('AAChange','.')
        SIFTprediction=INFOdict.get('SIFTprd','.')
        PP2prediction=INFOdict.get('PP2.hvar.prd','.')
        MAprediction=INFOdict.get('MutAprd','.')
        MTprediction=INFOdict.get('MutTprd','.')
        GERPscore=str(INFOdict.get('GERP','.'))
        CADDscore=str(INFOdict.get('CADDphred','.'))
        
        ## Get Sample Level Variant information
        QualityString=[ linelist[i].strip() for i in ColumnNumber ]
        QualityList=[ i.split(':') for i in QualityString ]
        GT=[ QualityList[i][0] for i in range(0,len(QualityString)) ]
        
        ## Check if all genotypes are present
        PassGeno=False
        if any(i not in BadGT  for i in GT):
            PassGeno=True
        
        PassGene=False
        if SearchGene in INFOstring:
            PassGene=True
        
        
        if PassGeno and PassGene:
            OutputList=linelist[0:5]+[GeneName,VariantFunction,VariantClass,AAchange,ExACFreq,KGFreq,ESPFreq,SIFTprediction,PP2prediction,MAprediction,MTprediction,GERPscore,CADDscore,SegDup]+GT+[FILTER,INFOstring]+QualityString
            OutputList= [ str(i) for i in OutputList ]
            OutputString="\t".join(OutputList)
            Output.write(OutputString+"\n")
            Outvcf.write(line)
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering Finished: "+TimeNow+"\n")
print "Filtering finished "+TimeNow
