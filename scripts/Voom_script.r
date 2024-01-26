
##############################################################
### Script for differential expression analysis using Voom ###
##############################################################

### Load limma-package (which includes Voom) ###
library(limma)
library(edgeR)     ## You also need this library to access the DGEList-function

### load data table ###
countTable = read.table("count-table.tsv",header=TRUE,row.names=1)   

### create condition ###
condition = factor( c("LNCaP","LNCaP","LNCaP","LNCaP","LNCaP","LNCaP",
			"RWPE","RWPE","RWPE","RWPE","RWPE","RWPE"))

### create design matrix ###
des = model.matrix(~-1+condition)
colnames(des) = levels(condition)

### define contrasts (which groups to compare) ###
cmat <- makeContrasts(LNCaP - RWPE, levels=des)

### Normalise count-table ###
dge <- DGEList(counts=countTable)
dge <- calcNormFactors(dge)

### Fit voom model ###
v <- voom(dge,design=des)
fit <- lmFit(v,design=des)
fit <- contrasts.fit(fit, cmat)
fit <- eBayes(fit)

### find differentially expressed transcripts ###
a <- decideTests(fit,adjust.method="fdr", p.value=0.05, lfc=0)

### summary of result-table ###
sma = summary(a)
dmm <- dim(countTable)
res <- topTable(fit,n=dmm[1],coef=1)

### write table with results from differential expression ###
write.table(res,file="Voom_diffexp.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

###########################################################################
