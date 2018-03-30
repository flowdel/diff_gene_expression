## Differential expression analysis 
library(Biobase)
library(BioNet)
library(GEOquery)
library(limma)
library(GEOquery)
library(genefilter)
library(impute)

# lymphoma data analysis
data(exprLym)
lymp <- exprLym

# make proper column names to match toptable 
fvarLabels(lymp) <- make.names(fvarLabels(lymp))

ex <- exprs(lymp)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(lymp) <- log2(ex) }

# set up the data and proceed with analysis
lymp$description <- lymp$Subgroup
design <- model.matrix(~ description + 0, lymp)
colnames(design) <- levels(lymp$Subgroup)
fit <- lmFit(lymp, design)
cont.matrix <- makeContrasts(ABC-GCB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# data table
tT_lym <- topTable(fit2, adjust="fdr", sort.by="logFC", number=3583)
tT_lym <- subset(tT_lym, select=c("adj.P.Val","logFC"))

# fit BUM sistribution
expressions <- impute.knn(exprs(exprLym))$data
t.test <- rowttests(expressions, fac = exprLym$Subgroup)
rownames(t.test) <- rownames(exprs(exprLym))
fb <- fitBumModel(t.test$p.value, plot = T)
hist(fb)
plot(fb)


# ESC-iPSC data analysis

# load series and platform data from GEO
gset <- getGEO("GSE71255", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "XXX000XXXXXXXX1111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")  
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# data table
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=45101)
tT <- subset(tT, select=c("adj.P.Val","logFC","Gene.symbol"))

# fit BUM distribution
fb <- fitBumModel(tT$adj.P.Val, plot = T)
hist(fb)
plot(fb)
