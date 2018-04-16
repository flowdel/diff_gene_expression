library(affy)
library(oligo)
library(pd.mogene.2.0.st)
library(genefilter)
library('mouse4302.db')
library('annotate')

setwd("/Volumes/Element/BI/project/GSE71255_RAW/")

data <- ReadAffy(filenames = c('GSM1831428_ESC1.CEL.gz', 'GSM1831429_ESC2.CEL.gz', 
                               'GSM1831440_OSKM2.CEL.gz', 'GSM1831439_OSKM1.CEL.gz'))
eset <- rma(data)
express <- exprs(eset)
eset$sample <- с(1, 1, 2, 2)
Index1 <-  which(eset$sample == 1)
Index2 <-  which(eset$sample == 2)
plot(rowMeans(express[,Index1]), rowMeans(express[,Index2]))
tt <- rowttests(express, factor(eset$sample))
design <- model.matrix(~factor(eset$sample))
fit <- lmFit(eset, design)
ebayes <- eBayes(fit)
tableTop <- topTable(ebayes,coef=2, number = 45101)

# annotation
genenames <- rownames(tableTop)
length(genenames)
annotation(eset)
geneID <- getEG(genenames, "mouse4302.db")
sym <- getSYMBOL(genenames, "mouse4302.db")
T1 <- data.frame(sym, tableTop)
T2 <- data.frame(geneID, T1)

write.table(x = express, file = 'expressions.csv', sep = '\t')
write.table(x = T2, file = 'table.csv', sep = '\t')
fb <- fitBumModel(T2$adj.P.Val)