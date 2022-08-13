specific_network <- function(sample_ID){
	sample_network_1 <- GRN[which(!is.na(match(GRN[,1],diff_mutation[,2]))),]
	sample_network <- sample_network_1[which(!is.na(match(sample_network_1[,2],rownames(resistant_exp)))),]
	sample_network[,5] <- rep("non-diff",nrow(sample_network))
	sample_network[which(!is.na(match(sample_network[,2],rownames(DEG)))),5] <- "diff"
	return(sample_network)
}