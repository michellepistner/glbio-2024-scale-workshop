library(Rcpp)
library(foreach)
library(doParallel)
library(MASS)
source("gsea_C.R")

### BASE GSEA FUNCTIONS ###

# Function that forms the basis of GSEA
# For GSEA-LFC scores is a vector of LFC
# inds are the indices in scores of the entities in the set of interest
gsea_base <- function(scores, inds) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	inds <- which(sorted_obj$ix%in%inds)
	inds <- sort(inds) 
	scores <- sorted_obj$x
	num_gt <- 0
	hit_scores <- abs(scores[inds])/sum(abs(scores[inds])) 
	return(gsea_score_C(hit_scores, inds, length(scores)) )
}

# Similar function but takes inds as a list
gsea_base_inds <- function(scores, inds, unweighted) {
	sorted_obj <- sort(scores, decreasing=T, index.return=T)
	scores <- sorted_obj$x
	gsea_scores <- rep(0, length(inds))
	for (i in 1:length(inds)) {
		idx <- sort(which(sorted_obj$ix%in%(inds[[i]])))
		hit_scores <- abs(scores[idx])/sum(abs(scores[idx]))
		if(unweighted) {
			hit_scores <- rep(1, length(hit_scores))
			hit_scores <- hit_scores /sum(hit_scores)
		}
		gsea_scores[i] <- gsea_score_C(hit_scores, idx, length(scores))
	}
	gsea_scores
}

### ENTITY LABEL PERMUTATIONS ###

gsea <- function(scores, inds, iters=5000, null_scores=NULL) {
	obs_score <- gsea_base(scores, inds)
	shuffled_scores <- c()
	if (!is.null(null_scores)) {
		shuffled_scores <- null_scores	
	} else {
		for (i in 1:iters) {
			shuffled_inds <- sample(1:length(scores), length(inds), replace=F)
			shuffled_scores <- c(shuffled_scores, gsea_base(scores, shuffled_inds))
		}
	}
	obs_sign <- sign(obs_score)
	shuffled_same_sign <- shuffled_scores[sign(shuffled_scores)==obs_sign]
	p_value <- sum(abs(shuffled_same_sign) >= abs(obs_score)) / length(shuffled_same_sign)
	return(list(obs_score=obs_score, p_value=p_value))
}

gsea_parallel <- function(scores, inds, iters=5000, cores=detectCores()) {
	cl <- makeCluster(cores)
	registerDoParallel(cl)
	clusterCall(cl, function() { source("gsea_C.R") })

	obs_scores <- gsea_base_inds(scores, inds, F)
	shuffled_scores <- foreach(i=1:iters, .combine=rbind, .export="gsea_base_inds") %dopar% {
		shuffled_inds <- lapply(inds, function(item) sample(1:length(scores), length(item), replace=F))
		gsea_base_inds(scores, shuffled_inds, F)
	}
	p_values <- c()
	for(i in 1:length(obs_scores)) {
		obs_sign <- sign(obs_scores[i])
		shuffled_same_sign <- (shuffled_scores[,i])[sign(shuffled_scores[,i])==obs_sign]
		p_values <- c(p_values, sum(abs(shuffled_same_sign) >= abs(obs_scores[i])) / length(shuffled_same_sign))
	}

	stopCluster(cl)
	gc()
	return(list(obs_scores=obs_scores, p_values=p_values))
}

lfc_sensitivity_testing <- function(scores, inds, iters, epsilons, cores=detectCores()) {
	cl <- makeCluster(cores)
	registerDoParallel(cl)

	clusterCall(cl, function() { source("gsea_C.R") })
	res <- foreach(j=1:length(inds), .export=c("gsea", "gsea_base"), .combine=c) %dopar% {
		for(e in sample(epsilons, replace=F)) {
			p_val <- gsea(scores+e, inds[[j]], iters=iters)$p_value
			if (is.na(p_val)) {
				next
			}
			if (p_val > 0.05) {
				return(F)
			}
		}
		return(T)
	}

	stopCluster(cl)
	gc()

	res
}

### FUNCTIONS FOR COLUMN LABEL PERMUTATIONS ###

ols_solution <- function(Y, X, lfc_column) {
	t(ginv(t(X)%*%X)%*%t(X)%*%t(Y))[,lfc_column]
}

draw_shuffled_col <- function(W, X, lfc_column) {
  X[,lfc_column] <- sample(X[,lfc_column], replace=F)
  ols_solution(W, X, lfc_column)
}

draw_perm_col <- function(W, X, X_shuf, lfc_column) {
  X[,lfc_column] <- X_shuf
  ols_solution(W, X, lfc_column)
}

score_mean_Euclidean_weights <- function(lfcs, all_lfcs) {
	me_weights <- c()
	for(lfc in lfcs) {
		me_weights <- c(me_weights, (1/(length(all_lfcs)-1)) * sum(abs(lfc-all_lfcs)))
	}
	me_weights
}

### COLUMN LABEL PERMUTATIONS ###

gsea_parallel_matrix <- function(S, X, inds, lfc_column, iterations=2000, permutation_matrix=NULL, cores=detectCores(), score_func=NULL, unweighted=F) {
	cl <- makeCluster(cores)
	registerDoParallel(cl)

	# Calculate observed
	obs_lfc <- ols_solution(S, X, lfc_column)
	if (!is.null(score_func)) {
		obs_lfc <- score_func(obs_lfc, obs_lfc)
	}
	obs_scores <- gsea_base_inds(obs_lfc, inds, unweighted)
	
	if(!is.null(permutation_matrix)) {
	  iterations <- ncol(permutation_matrix) 
	}

	clusterCall(cl, function() { source("gsea_C.R") })
	res <- foreach(j=1:iterations, .combine=rbind, .export=c("draw_shuffled_col", "ols_solution", "gsea_base_inds", "draw_perm_col"), .packages="MASS") %dopar% {
		if(is.null(permutation_matrix)) {
			shuffled_lfc <- draw_shuffled_col(S, X, lfc_column)
		} else {
			shuffled_lfc <- draw_perm_col(S, X, permutation_matrix[,j], lfc_column) 
		}
		if (!is.null(score_func)) {
			shuffled_lfc <- score_func(shuffled_lfc, shuffled_lfc)
		}
		shuffled_scores <- gsea_base_inds(shuffled_lfc, inds, unweighted)
		shuffled_scores
	}

	p_values <- c()
	for(i in 1:ncol(res)) {
		shuff_filt <- res[,i][sign(res[,i])==sign(obs_scores[i])]
		p_values <- c(p_values, sum(abs(shuff_filt) >= abs(obs_scores[i]))/length(shuff_filt))
	}

	stopCluster(cl)
	gc()

	return(list(obs_scores=obs_scores, p_values=p_values))
}

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

calculate_W <- function(Y, X, lfc_column, epsilon) {
	W_par <- t(t(Y)/colSums(Y))
	W_perp <- log(apply(W_par, 2, function(col) 1/gm_mean(col)))
	error <- epsilon * X[,lfc_column]
	W <- log(W_par)
	W <- sweep(W, 2, W_perp, "+")
	W <- sweep(W, 2, error, "+")
	W
}

gsea_s <- function(S, X, inds, iterations=2000, permutation_matrix=NULL, cores=detectCores(), score_func=NULL, unweighted=F) {
	# Calculate observed
	obs_lfc <- ols_solution(S, X, 2)
	if (!is.null(score_func)) {
		obs_lfc <- score_func(obs_lfc, obs_lfc)
	}
	enrichment_score <- gsea_base_inds(obs_lfc, inds, unweighted)
	
	if(!is.null(permutation_matrix)) {
	  iterations <- ncol(permutation_matrix) 
	}

        res <- c()
    for(j in 1:iterations) {
		if(is.null(permutation_matrix)) {
			shuffled_lfc <- draw_shuffled_col(S, X, 2)
		} else {
			shuffled_lfc <- draw_perm_col(S, X, permutation_matrix[,j], 2) 
		}
		if (!is.null(score_func)) {
			shuffled_lfc <- score_func(shuffled_lfc, shuffled_lfc)
		}
		shuffled_scores <- gsea_base_inds(shuffled_lfc, inds, unweighted)
		res <- rbind(res, shuffled_scores)
	}

	p_value <- c()
	for(i in 1:ncol(res)) {
		shuff_filt <- res[,i][sign(res[,i])==sign(enrichment_score[i])]
		p_value <- c(p_value, sum(abs(shuff_filt) >= abs(enrichment_score[i]))/length(shuff_filt))
	}

    ret_m <- cbind(enrichment_score, p_value)
    row.names(ret_m) <- names(inds)
	return(ret_m)
}

