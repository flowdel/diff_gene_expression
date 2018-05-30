library(BioNet)
library(DLBCL)
library(ggplot2)
library(Hmisc)
library(sqldf)
library(igraph)
library(dplyr)
library(tidyr)
library(MASS)
library(limma)

# load and clean data 
data(exprLym)
data(interactome)

data <- exprs(exprLym)
colnames(data) <- exprLym$Subgroup
data[which(is.na(data))] <- apply(data[c(row(data)[c(which(is.na(data)))]),],1,function(x) x[c(which(is.na(x)))] <- mean(x, na.rm = TRUE))

# separate two lymphomas
gcb <- subset(data,select = colnames(data)=='GCB')
abc <- subset(data,select = colnames(data)=='ABC')

# create matrix and find correlations
gcb_m <- matrix(data=t(gcb),nrow=112,ncol=3583)
abc_m <- matrix(data=t(abc),nrow=82,ncol=3583)
colnames(gcb_m) <- rownames(gcb)
colnames(abc_m) <- rownames(abc)
cor_data_gcb <- rcorr(gcb_m)
cor_data_abc <- rcorr(abc_m)

# eject significant gene correlations
rcorr_p_gcb <- which(cor_data_gcb[["P"]]<0.05 & cor_data_gcb[["r"]]>0.5)
rcorr_p_abc <- which(cor_data_abc[["P"]]<0.05 & cor_data_abc[["r"]]>0.5)

col_rcorr_p_gcb <- rcorr_p_gcb%%3583
row_rcorr_p_gcb <- rcorr_p_gcb%/%3583+1

col_rcorr_p_abc <- rcorr_p_abc%%3583
row_rcorr_p_abc <- rcorr_p_abc%/%3583+1

row_rcorr_p_gcb[which(col_rcorr_p_gcb==0)] <- row_rcorr_p_gcb[which(col_rcorr_p_gcb==0)] - 1
col_rcorr_p_gcb[which(col_rcorr_p_gcb == 0)] <- c(3583)
rcorr_p_gcb_frame <- as.data.frame(c(rownames(cor_data_gcb$P)[col_rcorr_p_gcb]),c(colnames(cor_data_gcb$P)[row_rcorr_p_gcb]))

row_rcorr_p_abc[which(col_rcorr_p_abc==0)] <- row_rcorr_p_abc[which(col_rcorr_p_abc==0)] - 1
col_rcorr_p_abc[which(col_rcorr_p_abc == 0)] <- c(3583)
rcorr_p_abc_frame <- as.data.frame(c(rownames(cor_data_abc$P)[col_rcorr_p_abc]),c(colnames(cor_data_abc$P)[row_rcorr_p_abc]))


# find correlated genes that have edge in interactome graph
g <- igraph.from.graphNEL(interactome, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
compg.edges <- as.data.frame(get.edgelist(g))

a1 <- rcorr_p_gcb_frame
a2 <- compg.edges
a1Ina2_gcb <- sqldf('SELECT * FROM a1 INTERSECT SELECT * FROM a2')
head(a1Ina2_gcb)

b1 <- rcorr_p_abc_frame
b2 <- compg.edges
a1Ina2 <- sqldf('SELECT * FROM b1 INTERSECT SELECT * FROM b2')
head(b1Inb2_abc)

# save interactions in different lymphomas if it's necessary
m1 <- as.matrix(as.data.frame(cor_data_gcb$r))
df_gcb <- data.frame(row=rownames(m1)[row(m1)], col=colnames(m1)[col(m1)], corr=c(m1))
colnames(df_gcb) <- c('V1','V2','corr')
df_gcb <- df_gcb[!duplicated(df_gcb), ]

m2 <- as.matrix(as.data.frame(cor_data_abc$r))
df_abc <- data.frame(row=rownames(m2)[row(m2)], col=colnames(m2)[col(m2)], corr=c(m2))
colnames(df_abc) <- c('V1','V2','corr')
df_abc <- df_abc[!duplicated(df_abc), ]

d_gcb <- right_join(df_gcb, a1Ina2_gcb, by = c("V1",'V2'))
colnames(d_gcb) <- c('from','to','corr')

d_abc <- right_join(df_abc, b1Inb2_abc, by = c("V1",'V2'))
colnames(d_abc) <- c('from','to','corr')

write.table(d_gcb, file='corr_gcb.csv',sep=';')
write.table(d_abc, file='corr_abc.csv',sep=';')


# try to find non-linear correlations with chi_squared function
chi_sq <- function(t1, t2){
  x <- sort(t1)
  y <- sort(t2)
  n <- 3
  chunk <- function(x, n) split(x, sort(rank(x) %% n))
  c1 <- chunk(x,n)
  c2 <- chunk(y,n)
  mat <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:length(t1)){
    if(t1[i] %in% c1$`0` & t2[i] %in% c2$`0`){
      mat[1,1]<-mat[1,1]+1
    }
    if(t1[i] %in% c1$`1` & t2[i] %in% c2$`0`){
      mat[1,2]<-mat[1,2]+1
    }
    if(t1[i] %in% c1$`2` & t2[i] %in% c2$`0`){
      mat[1,3]<-mat[1,3]+1
    }
    if(t1[i] %in% c1$`0` & t2[i] %in% c2$`1`){
      mat[2,1]<-mat[2,1]+1
    }
    if(t1[i] %in% c1$`1` & t2[i] %in% c2$`1`){
      mat[2,2]<-mat[2,2]+1
    }
    if(t1[i] %in% c1$`2` & t2[i] %in% c2$`1`){
      mat[2,3]<-mat[2,3]+1
    }
    if(t1[i] %in% c1$`0` & t2[i] %in% c2$`2`){
      mat[3,1]<-mat[3,1]+1
    }
    if(t1[i] %in% c1$`1` & t2[i] %in% c2$`2`){
      mat[3,2]<-mat[3,2]+1
    }
    if(t1[i] %in% c1$`2` & t2[i] %in% c2$`2`){
      mat[3,3]<-mat[3,3]+1
    }
  }
  mat1 <- as.table(mat)
  return(chisq.test(mat1))
}


# simulate data with normal distribution to use it for finding significant chi_sq dependencies
sim_data <- matrix(0,nrow=3583,ncol=112)
for (i in 1:nrow(sim_data)){
  sim_data[i,] <- rnorm(112,0,1)
}

# start chi_sq function
mat_gcb_p <- matrix(-1,nrow=3583,ncol=3583)
mat_gcb_chi <- matrix(-1,nrow=3583,ncol=3583)
for (i in 1:(nrow(gcb)-1)){
  for (j in (i+1):nrow(gcb)){
    val <- chi_sq(gcb[i,], gcb[j,])
    mat_gcb_chi[i,j] <- val$statistic
    mat_gcb_p[i,j] <- val$p.value
  }
}

# save chi_sq function results
write.table(mat_gcb_chi, file='chi_gcb.csv',sep=';')
write.table(mat_gcb_p, file='chi_p_gcb.csv',sep=';')

# comparison of chi_data and rcorr_data for GCB lymphoma
data_chi <- read.csv('chi_p_gcb.csv', sep=';')
rownames(data_chi) <- rownames(gcb)
colnames(data_chi) <- rownames(gcb)

chi_val <- which(data_chi<0.00000004 & data_chi>-1)
col_chi <- chi_val%%3583
row_chi <- chi_val%/%3583+1
row_chi[which(col_chi==0)] <- row_chi[which(col_chi==0)] - 1
col_chi[which(col_chi == 0)] <- c(3583)
chi_frame <- as.data.frame(c(rownames(data_chi)[col_chi]),c(colnames(data_chi)[row_chi]))

df_chi <- as.data.frame(matrix(nrow=1869, ncol=2))
df_chi$V1 <- rownames(chi_frame)
df_chi$V2 <- as.character(chi_frame$`c(rownames(data_chi)[col_chi])`)
df_chi_2 <- as.data.frame(matrix(ncol=2, nrow=3738))
df_chi_2[1:1869,1] <- df_chi$V1
df_chi_2[1:1869,2] <- df_chi$V2
df_chi_2[1870:3738,1] <- df_chi$V2
df_chi_2[1870:3738,2] <- df_chi$V1

df_corr <- as.data.frame(matrix(nrow=33312, ncol=2))
df_corr$V1 <- rownames(rcorr_p_gcb_frame)
df_corr$V2 <- as.character(rcorr_p_gcb_frame$`c(rownames(cor_data_gcb$P)[col_rcorr_p_gcb])`)

a1 <- df_chi_2
a2 <- df_corr

chisq_In_corr <- sqldf('SELECT * FROM a1 INTERSECT SELECT * FROM a2')
chisq_Not_In_corr <- sqldf('SELECT * FROM a1 EXCEPT SELECT * FROM a2')


# data simulation with covariation matrix
cov_gcb <- cov(gcb_m)
cov_abc <- cov(abc_m)
# lists of mean values
z1 <- rep(0, 3583)
z2 <- rep(0.2, 3583)

sim_data_1 <- mvrnorm(112, mu = z1, Sigma = cov_gcb)
sim_data_2 <- mvrnorm(82, mu = z1, Sigma = cov_abc)
rownames(sim_data_1) <- rep('GCB', 112)
rownames(sim_data_2) <- rep('ABC', 82)
sim_1 <- t(sim_data_1)
sim_2 <- t(sim_data_2)

# convert data to expression set
d <- cbind(sim_1, sim_2)
colnames(d) <- as.character(c(1:194))
fac <- rep('GCB', 112)
fac <- c(fac, rep('ABC', 82))
fac <- as.factor(fac)
d1 <- t(d) 
d1 <- ExpressionSet(d)
d1$description <- fac

# analyse expression set
design <- model.matrix(~ description + 0, d1)
colnames(design) <- levels(fac)
fit <- lmFit(d1, design)
cont.matrix <- makeContrasts(GCB-ABC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=3583)
fb <- fitBumModel(tT$P.Value, plot = T)



