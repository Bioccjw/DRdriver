differential_exp <- function(sample_ID){
	library(edgeR)
	trans_ID <- sample_transformation(sample_ID)
	sample_resistant_exp <- as.matrix(resistant_exp[,match(trans_ID,colnames(resistant_exp))])
	rownames(sample_resistant_exp) <- rownames(resistant_exp)
	colnames(sample_resistant_exp) <- trans_ID
	merge_exp <- cbind(sample_resistant_exp,sensitive_exp)
	merge_exp <- as.matrix(merge_exp)
	group <- c(1,rep(2,ncol(merge_exp)-1))
	y <- DGEList(merge_exp,group)
	keep <- rowSums(cpm(y)>1) >= 2 
	y <- y[keep, , keep.lib.sizes=FALSE]
	y <- calcNormFactors(y)
	design <- model.matrix(~group)
	y <- estimateDisp(y,design)
	filt <- glmQLFit(y,design)
	lrt <- glmQLFTest(filt,coef=2)
	gene_condition <- topTags(lrt,n=60000)
	result <- gene_condition$table
	DEG <- result[which(result[,4]<0.05),]
    return(DEG)
}