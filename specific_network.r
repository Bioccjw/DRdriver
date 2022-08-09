specific_network <- function(sample_ID){
    diff_mutation <- differential_mutation(sample_ID)
	DEG <- differential_exp(sample_ID)
	sample_network_1 <- GRN[which(!is.na(match(GRN[,1],diff_mutation[,2]))),]
	sample_network <- sample_network_1[which(!is.na(match(sample_network_1[,2],rownames(resistant_exp)))),]
	sample_network[,3] <- rep("non-diff",nrow(sample_network))
	sample_network[which(!is.na(match(sample_network[,2],rownames(DEG)))),3] <- "diff"
	return(sample_network)
}
