library(affy)
library(oligo)
library(pd.mogene.2.0.st)
library(genefilter)
library(annotate)
library(annaffy)

# annotation database
library('mouse4302.db')
# for collapseBy
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")

# first you must download GSE data from GEO
# then load GSE data
setwd("/GSE71255_RAW/")
data <- ReadAffy(filenames = c('GSM1831426_MEF2.CEL.gz', 'GSM1831427_MEF3.CEL.gz', 
                               'GSM1831428_ESC1.CEL.gz', 'GSM1831429_ESC2.CEL.gz',
                               'GSM1831430_ESC3.CEL.gz', 'GSM1831440_OSKM2.CEL.gz', 
                               'GSM1831439_OSKM1.CEL.gz', 'GSM1831441_OSKM3.CEL.gz',
                               'GSM1831442_OSKM4.CEL.gz'))

# convert data into expression set
eset <- affy::rma(data)

# define annotation database (then download it)
annotation(eset)

# annotate expression set
probenames <- rownames(eset)
fData(eset)$gene_id <- getEG(probenames, "mouse4302.db")
fData(eset)$gene_symbol <- getSYMBOL(probenames, "mouse4302.db")

# filter data
es <- collapseBy(eset, fData(eset)$gene_id, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]
es_top <- es[1:12000, ]
es_top$sample <- c("MEF", "MEF", "ESC", "ESC", "ESC", "OKSM", "OKSM", "OKSM", "OKSM")

# if it's necessary, normalize data
#exprs(es_top) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")






