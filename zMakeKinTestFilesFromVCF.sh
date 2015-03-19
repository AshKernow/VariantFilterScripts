#!/bin/bash

VcfFil=March06_2015_ClinicalRoR.rawvariants.Ann.hardfiltered.vcf
OutNam=${VcfFil/.vcf/}

qsss -c "vcftools --vcf $VcfFil --relatedness2 --out $OutNam"
qsss -c "vcftools --vcf $VcfFil --plink --out $OutNam"
plink --file $OutNam --make-bed --split-x hg19 --out $OutNam
plink --bfile $OutNam --impute-sex --make-bed --out $OutNam.sexcheck
plink --bfile $OutNam --genome --out $OutNam
mv $OutNam.sexcheck.sexcheck $OutNam.sexcheck
rm -rf $OutNam.ped $OutNam.map *nosex *.sexcheck.*


