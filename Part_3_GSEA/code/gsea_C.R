library(Rcpp)

# C++ optimized function for computing GSEA score.
# It is assumed that the hit scores / inds are in sorted order
# hit_scores is a vector of the fraction weights of each gene e.g., c(0.5,0.3,0.2)
# num_genes is the total number of genes
cppFunction('
double gsea_score_C(NumericVector hit_scores, NumericVector inds, int num_genes) {
	double score = 0;
	double max_score = 0;
	double prev = 0;
	double miss_score = -1.0/(num_genes-inds.size());
	for(int i=0;i<inds.size();++i) {
		score += (inds[i]-prev-1)*miss_score;
		if (fabs(score) > fabs(max_score)) {
			max_score = score;
		}

		score += hit_scores[i];
		prev = inds[i];
		if (fabs(score) > fabs(max_score)) {
			max_score = score;
		}
	}
	return max_score;
}
')

