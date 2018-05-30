library(igraph)
library(BioNet)
library(knitr)
library(mwcsr)
library(mcmcRanking)
library(rlist)

# first you must download InWeb PPI graph
# then load InWeb PPI graph
int <- as.data.frame(read.csv('/inwebIM_ppi.txt', sep='\t'))
g <- graph_from_data_frame(int, directed = FALSE)

#ESC vs OKSM
# analyse obtained in "CEL_to_expression_set.R" data and get p-values
design <- model.matrix(~0+sample, data=pData(es_top))
fit <- lmFit(es_top, design)
fit2 <- contrasts.fit(fit, makeContrasts(sampleESC-sampleOKSM, 
                                         levels=design))
ebayes <- eBayes(fit2)
tableTop_ESCvsOKSM <- topTable(ebayes,adjust.method="BH",number = Inf)
pval_ESCvsOKSM <- as.vector(tableTop_ESCvsOKSM$P.Value)
names(pval_ESCvsOKSM) <- toupper(tableTop_ESCvsOKSM$gene_symbol)

# remove some genes
pval_ESCvsOKSM <- list.remove(pval_ESCvsOKSM, 'UBC')
pval_ESCvsOKSM <- list.remove(pval_ESCvsOKSM, 'SUMO1')
pval_ESCvsOKSM <- list.remove(pval_ESCvsOKSM, 'SUMO2')
pval_ESCvsOKSM <- list.remove(pval_ESCvsOKSM, 'SUMO3')
pval_ESCvsOKSM <- list.remove(pval_ESCvsOKSM, 'APP')
length(pval_ESCvsOKSM)

logFC_ESCvsOKSM <- as.vector(tableTop_ESCvsOKSM$logFC)
names(logFC_ESCvsOKSM) <- toupper(tableTop_ESCvsOKSM$gene_symbol)
logFC_ESCvsOKSM <- list.remove(logFC_ESCvsOKSM, 'UBC')

# get subgraph
graph_ESCvsOKSM <- largestComp(induced.subgraph(g, which(V(g)$name %in% names(pval_ESCvsOKSM))))

# check if it's connected subgraph
is.connected(graph_MEFvsOKSM) 

# rmwcs analysis
solve <- rmwcs(timelimit = 30, max_iterations = 30)
fb <- fitBumModel(pval_ESCvsOKSM, plot = F)
# you can change fdr
scores <- scoreNodes(graph_ESCvsOKSM, fb, fdr = 0.00000001)
V(graph_ESCvsOKSM)$score <- scores
# find module
m <- solve(graph_ESCvsOKSM)
# visualise module
plotModule(m, scores = scores, diff.expr = logFC_ESCvsOKSM)


# mcmcRanking analysis
# set p-values
pval_from_graph_ESCvsOKSM <- pval_ESCvsOKSM[which(names(pval_ESCvsOKSM) %in% names(V(graph_ESCvsOKSM)))]
pval_from_graph_ESCvsOKSM <- pval_from_graph_ESCvsOKSM[order(match(names(pval_from_graph_ESCvsOKSM),names(V(graph_ESCvsOKSM))))]
V(graph_ESCvsOKSM)$pval <- pval_from_graph_ESCvsOKSM
# set likelihood
V(graph_ESCvsOKSM)$likelihood <- exp(scores1)
# find module
# you can change module size
x <- mcmc_sample(graph = graph_ESCvsOKSM, module_size = 69, times = 100, iter = 2e5)
mod <- induced.subgraph(graph_ESCvsOKSM, x$name[which(x$mat[1,])])
# visualise module
test.layout <- layout_(mod,with_dh(weight.edge.lengths = 1/1000))
plot(mod, vertex.color='yellow',layout = test.layout)

