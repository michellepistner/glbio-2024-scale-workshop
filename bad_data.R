#### Code and data for 'bad data' figures
library(DESeq2)
library(ALDEx2)


# Dataset 1, the synthetic dataset

Y <- t(rdat[,-1])
Y[1:5,1:5]
colnames(Y) <- rownames(rdat)
colnames(Y)[1:100] <- paste('pre',colnames(Y)[1:100], sep="")
colnames(Y)[101:200] <- paste('post',colnames(Y)[101:200], sep="")

conds <- as.character(rdat[,1])
conds

colData <- cbind(conds, rep(1,length(conds)))
colnames(colData) <- c('conditions','batch')
rownames(colData) <- colnames(Y)

# analyze with DESeq2
dds <- DESeqDataSetFromMatrix(countData = Y,
                              colData = colData,
                              design= ~conditions)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="conditions_Pre_vs_Post")

# display the adjusted p values
# note that almost all of them are significant
hist(res$padj, xlim=c(0,1), breaks=19)
abline(v=0.05, col='red', lty=2)

#generate an MA plot
plotMA(res)

sig <- res$padj < 0.05

# generate an effect plot
disp <- dds@rowRanges@elementMetadata@listData$SE_conditions_Pre_vs_Post*sqrt(20)
diff <- -res$log2FoldChange

plot(disp,diff, xlab="SD change", ylab="log2FC", pch=19, cex=0.5, col='grey')
points(disp[sig],diff[sig],pch=0, cex=0.8)
abline(h=0)