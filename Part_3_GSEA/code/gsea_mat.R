library(fgsea)
source("gsea.R")

pathways <- gmtPathways("c2.all.v7.4.symbols.gmt")
Y <- as.matrix(read.table("breast_tumor_counts.txt", row.names=1))
metadata <- read.table("breast_tumor_metadata.txt", header=T)
metadata[,2] <- factor(metadata[,2], level=c("Healthy", "Tumor"))
head(metadata)
row.names(metadata) <- NULL
X <- model.matrix(~condition, metadata)
path_inds <- list()
path_names <- list()
for (p in names(pathways)) {
	if(!grepl("KEGG_", p)) {
		next
	}
	inds <- match(pathways[[p]], row.names(Y))
        inds <- inds[!is.na(inds)]
        if((length(inds)<15)|(length(inds)>500)) {
                next
        }
        path_inds[p] <- list(inds)
        path_names[[p]] <- pathways[[p]][pathways[[p]]%in%row.names(Y)]
}

gsea_s.error <- function(W, X, pathways, iterations=500, epsilon=c(0, 0.5)) {
    all_gsea_res <- c()
    for(epsilon_perp in epsilon) {
            noise <-c(rep(0, 92), rep(epsilon_perp, 92))
            noise_adj_W <- sweep(W, 2, noise, "+")
            gsea_res <- gsea_s(noise_adj_W, X, pathways, iterations=iterations)
            gsea_res <- cbind(epsilon_perp, gsea_res)
            all_gsea_res <- rbind(all_gsea_res, gsea_res)
    }
    return(all_gsea_res)
}

W <- calculate_W(Y, X, 2, 0)
gsea_res <- gsea_s(W, X, path_inds, iterations=250)
gsea_lfcs_sens_res <- gsea_s.error(W, X, path_inds, iterations=250, epsilon=c(0,0.5))
head(gsea_lfcs_sens_res)

