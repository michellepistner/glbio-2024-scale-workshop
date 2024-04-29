## Analysis of Vandputte et al (2017) "Quantitative Microbiome Profiling..." Nature.
## Helper function for computing the LFC----------------------------------------
aldex.lfc <- function(clr){
  # Use clr conditions slot instead of input
  if (is.vector(clr@conds)) {
    conditions <- clr@conds
  } else if (is.factor(clr@conds)) {
    if (length(levels(clr@conds) == 2)) {
      conditions <- clr@conds
    }
  } else if (is.matrix(clr@conds)){
    stop("currently does not support > 2 conditions.")
  } else {
    stop("please check that the conditions parameter for aldex.clr is correct.")
  }
  
  nr <- numFeatures(clr) # number of features
  rn <- getFeatureNames(clr) # feature names
  mc.s <- numMCInstances(clr)
  p <- length(conditions)
  # ---------------------------------------------------------------------
  
  # sanity check to ensure only two conditons passed to this function
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  sets <- levels
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  
  ## we need to fill an array and calculate these statistics over the array
  W.est <- array(NA, dim = c(nr, p, mc.s))
  
  for(i in 1:p){
    W.est[,i,] <- getMonteCarloReplicate(clr, i)
  }
  print(paste0("Condition A is ", sets[1], " and condition B is ", sets[2], "."))
  lfc <- rep(NA, mc.s)
  lfc <- apply(W.est, MARGIN = 3, FUN = function(mat, setA, setB){rowMeans(mat[,setA]) - rowMeans(mat[,setB])}, setA =setA, setB = setB)
  return(data.frame("lfc" = rowMeans(lfc), "sd" = apply(lfc,1,sd), "p2.5" = apply(lfc,1,quantile, probs = c(0.025)), "p97.5" = apply(lfc,1,quantile, probs = c(.975))))
}

#-------------------------------------------------------------------------------

# Analysis----------------------------------------------------------------------

## Loading libraries
library(phyloseq) #
library(tidyverse) #
library(ALDEx2) #
library(gghighlight)
library(cowplot)
library(ggplot2)
library(ggpattern)
library(magrittr)
library(stringi)
library(phyloseq)
library(data.table)
library(DESeq2)
library(edgeR)
library(latex2exp)
library(ggrepel)
set.seed(2024)

## Reading in the phyloseq object
phylo <- readRDS(file.path("Part 2 - Updates to ALDEx2", "data","vandputte", "vandputte_phyloseq.RDS"))

## Investigating the sample data
head(sample_data(phylo))

## Plotting the flow cytometry data
sample_data(phylo) %>% 
  data.frame() %>%
  dplyr::select(Health.status, CellCount) %>%
  ggplot(aes(x = Health.status, y = CellCount)) +
  geom_boxplot() +
  xlab("Cohort") +
  ylab("Flow Cytometry Measurement")


## Plotting the implied scale model for the original ALDEx2 model
clr <- aldex.clr(t(otu_table(phylo)), sample_data(phylo)$Health.status, mc.samples = 1000, denom = "all", gamma = 1e-3)
clr.dat <- data.frame("group" = sample_data(phylo)$Health.status, "scale" = clr@scaleSamps[,1])
ggplot(clr.dat, aes(x=group, y=scale)) + 
  geom_boxplot() +
  xlab("Cohort") +
  ylab("Implied Value of Scale")


#-------------------------------------------------------------------------------

## ALDEx2 models----------------------------------------------------------------
## First aldex2
Y <- t(otu_table(phylo)) 
X <- sample_data(phylo)$Health.status

##Default scale model
aldex_fit <- aldex(Y,X,mc.samples = 1000, gamma = .5)

tax_default = aldex_fit %>% 
  rownames_to_column("category") %>%
  dplyr::select(category, effect, we.ep, we.eBH) %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)

##Gold standard flow cytometry model

## Choosing scale based off of the estimated technical variation from the other methods. See "Data" for the calculation in the excel file "41586_2017_BFnature24460_MOESM11_ESM.xlsx"
## Doing a sensitivity analysis as well
sen_res <- list()
scale_mean <- log2(sample_data(phylo)$CellCount)
gamma <-  c(1e-3,.1,.25,.5,0.7,1,1.5,2, 3,4,5,6,7,8,9,10)
for(j in 1:length(gamma)){
  scale_var <- rep(gamma[j], 95)
  
  scale_samples <- matrix(NA, nrow = 95, ncol = 1000)
  for(i in 1:95){
    scale_samples[i,] <- 2^rnorm(1000, scale_mean[i], scale_var[i])
  }
  clr <- aldex.clr(Y,X, gamma = scale_samples,mc.samples = 1000)
  res <- aldex.ttest(clr)
  res <- cbind(res, aldex.lfc(clr))
  sen_res[[j]] <- res
}

## Creating the graph data frame
effect <- c()
effect.low <- c()
effect.high <- c()
lfc.sd <- c()
gam_used <- c()
we.eBH <- c()
tax <- c()

for(i in 1:length(sen_res)){
  gam_used <- c(gam_used, rep(gamma[i], nrow(sen_res[[i]])))
  effect <- c(effect, sen_res[[i]]$lfc)
  effect.low <- c(effect.low, sen_res[[i]]$p2.5)
  effect.high <- c(effect.high, sen_res[[i]]$p97.5)
  lfc.sd <- c(lfc.sd, sen_res[[i]]$sd)
  we.eBH <- c(we.eBH, sen_res[[i]]$we.eBH)
  tax <- c(tax, rownames(sen_res[[i]]))
}

graph.df <- data.frame("gamma" = gam_used, "effect" = effect, "we.eBH" = we.eBH, "Sequence" = tax, "low" = effect.low, "high" = effect.high, "sd" = lfc.sd)

tax_merge <- as.data.frame(tax_table(phylo)) %>%
  rownames_to_column("Sequence") %>%
  dplyr::select(Sequence, Genus)

graph.df <- graph.df %>%
  mutate(effect = -1*effect/sd) %>%
  plyr::join(tax_merge, by = "Sequence") %>%
  mutate(Genus = ifelse(Genus == "Escherichia/Shigella", "Escherichia", Genus)) %>%
  mutate(label = ifelse((we.eBH < 0.05), Genus, NA)) %>%
  group_by(label) %>%
  mutate(gamma.max = max(gamma)) %>%
  mutate(label = ifelse((!is.na(label)) & (gamma == gamma.max), label, NA)) %>%
  mutate(label = ifelse((!is.na(label)) & (gamma <= 6), label, NA)) %>%
  mutate(label = ifelse(label %in% c("Anaerobutyricum", "Sutterella", "Parabacteroides"), label, NA))

gam_diag <- ggplot(graph.df, aes(x = gamma, y = effect, group = Genus)) +
  geom_line () +
  gghighlight(we.eBH < 0.05) +
  ylab("Standardized Log Fold Change") +
  xlab(expression(gamma)) +
  theme(text = element_text(size=28)) +
  theme_bw() +
  geom_label_repel(aes(label = label),
                   na.rm = TRUE, size=3, direction='y') 


## Choosing scale based off of the estimated technical variation from the other methods. See "Data" for the calculation in the excel file "41586_2017_BFnature24460_MOESM11_ESM.xlsx"
tax_scale = sen_res[[5]] %>% 
  rownames_to_column("category") %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)


## Informed model
scale.cd <- 2^matrix(rnorm(1000*29, mean = log2(.7), sd = .125), nrow = 29)
scale.control <- 2^matrix(rnorm(1000*66, mean = log2(1), sd = .125), nrow = 66)

scale.inf <- rbind(scale.cd, scale.control)
aldex_informed <- aldex(Y,X,mc.samples = 1000, gamma = scale.inf)

tax_informed = aldex_informed %>% 
  rownames_to_column("category") %>%
  dplyr::select(category, effect, we.ep, we.eBH) %>%
  mutate(pval = we.ep) %>%
  mutate(padj = we.eBH) %>%
  dplyr::filter(padj < 0.05)

#-------------------------------------------------------------------------------

## Other methods ---------------------------------------------------------------

## DeSeq2
coldata <- matrix(sample_data(phylo)$Health.status, ncol = 1)
colnames(coldata) <- "Condition"
dds <- DESeqDataSetFromMatrix(countData=as.matrix(Y@.Data + 1),
                              colData=coldata,
                              design = ~Condition)
dds <- DESeq(dds)
res <- results(dds)
tax_deseq <- res %>% 
  as.data.frame() %>% 
  rownames_to_column("category") %>% 
  dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
  mutate(sig = ifelse(padj < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)


##EdgeR

Y <- t(otu_table(phylo))
X <- sample_data(phylo)$Health.status
conds <- sample_data(phylo)$Health.status 
coldata <- ifelse(conds == "CD", 0, 1)
coldata <- matrix(coldata, ncol = 1)

y <- DGEList(counts=as.matrix(Y@.Data +1),group=c(coldata[,1]))
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tax_edgeR <- topTags(et, n=1600, p.value = 1) %>%
  as.data.frame() %>%
  rownames_to_column("category") %>%
  mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

#-------------------------------------------------------------------------------

## vandpuette's qmp method------------------------------------------------------

# function from https://github.com/raeslab/QMP/blob/master/QMP.R
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[row.names(cnv_corrected_abundance_table),]) # make sure the order of the samples is the same in both files    
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}
cellCounts <- data.frame(sample_data(phylo)) %>%
  dplyr::select(CellCount)
QMP_df <- rarefy_even_sampling_depth(otu_table(phylo),cellCounts)

## running a wilcoxon test on the qmp data
p.qmp <- rep(NA, ncol(QMP_df))
lfc.qmp <- rep(NA, ncol(QMP_df))
conds <- sample_data(phylo)$Health.status
for(i in 1:length(p.qmp)){
  x <- log2(QMP_df[which(conds == "CD"),i]+0.5)
  y <- log2(QMP_df[which(conds == "Control"),i]+0.5)
  lfc.qmp[i] <- mean(x) - mean(y)
  tmp.test <- t.test(x,y)
  p.qmp[i] <- tmp.test$p.value
}
padj.qmp <- p.adjust(p.qmp, method = "BH")

tax_QMP <- data.frame("category" = colnames(otu_table(phylo)),
                      "pval" = p.qmp,
                      "padj" = padj.qmp) %>%
  mutate(sig = ifelse(padj.qmp <= 0.05, TRUE, FALSE)) %>%
  dplyr::filter(sig == TRUE)

#-------------------------------------------------------------------------------

## Plotting---------------------------------------------------------------------

sig.values = c(tax_edgeR$category,tax_deseq$category, tax_scale$category, tax_informed$category, tax_QMP$category, tax_default$category) %>% unique

###Some light processing to make it more useful

truth.pos = tax_scale$category

##Generating the grid plot
q=length(sig.values)

sig.df = data.frame("Sequence" = rep(sig.values,6))
sig.df = sig.df %>%
  mutate(true.pos = ifelse(Sequence %in% truth.pos, 1, 0)) %>%
  mutate(Model = c( rep("ALDEx2 (Gold Standard)", q), rep("ALDEx2 (Informed)", q), rep("ALDEx2 (Default)", q), rep("QMP", q),rep("DESeq2", q), rep("edgeR", q))) %>%
  mutate(sigcode = c( ifelse(sig.values %in% tax_scale$category, 1, 0),ifelse( sig.values %in% tax_informed$category, 1, 0),ifelse( sig.values %in% tax_default$category, 1, 0), ifelse( sig.values %in% tax_QMP$category, 1, 0),ifelse( sig.values %in% tax_deseq$category, 1, 0), ifelse( sig.values %in% tax_edgeR$category, 1, 0))) %>%
  mutate(res = ifelse(true.pos == 1 & sigcode == 1, "TP", NA)) %>%
  mutate(res = ifelse(true.pos == 0 & sigcode == 0, "TN", res)) %>%
  mutate(res = ifelse(true.pos == 1 & sigcode == 0, "FN", res)) %>%
  mutate(res = ifelse(true.pos == 0 & sigcode == 1, "FP", res)) %>%
  mutate(sigcode = factor(sigcode, levels = list("Non. Sig."="0", "Sig."="1")))


sig.df$Model = factor(sig.df$Model, levels=c( "edgeR", "DESeq2","QMP","ALDEx2 (Default)", "ALDEx2 (Informed)", "ALDEx2 (Gold Standard)"))

tax_merge <- as.data.frame(tax_table(phylo)) %>%
  rownames_to_column("Sequence") %>%
  dplyr::select(Sequence, Genus)

sig.df <- sig.df %>%
  plyr::join(tax_merge, by = "Sequence") %>%
  mutate(Genus = ifelse(Genus == "Escherichia/Shigella", "Escherichia", Genus)) 
##No color labels
p1 = ggplot(sig.df, aes(x=Genus, y=Model)) +
  geom_tile_pattern(aes(fill=res, pattern = res), color="black",pattern_fill = 'black', pattern_colour  = 'black', pattern_density = 0.015) +
  theme_minimal(18) +
  labs(title = "") +
  theme(panel.grid = element_blank(), 
        legend.title=element_blank(),
        legend.key.size = unit(.825, "cm"),
        text = element_text(size=20)) +
  scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "stripe", FN = "stripe")) +
  scale_fill_manual(values= c("white", "lightgrey", "white", "grey")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1
#-------------------------------------------------------------------------------
