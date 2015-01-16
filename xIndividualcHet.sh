#!/bin/bash
#$ -l mem=4G,time=1:: -cwd -N IndcHet

R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

fils <- list.files(pattern="IcHet.*.tsv")

for(i in 1:length(fils)){
    dat <- read.delim(fils[i])
    dups <- unique(dat[duplicated(dat[,"Gene"]),"Gene"])
    dat <- dat[dat[,"Gene"]%in%dups,]
    write.table(dat, fils[i], quote=F, row.names=F, sep="\t")
}
RSCRIPT
