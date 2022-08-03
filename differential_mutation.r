differential_mutation <- function(sample_ID){
    sample_resistance_mutation <- resistant_mutation[which(resistant_mutation[,1]==sample_ID),]
	specific_mutation <- setdiff(sample_resistance_mutation[,3],sensitive_mutation[,3])
	resistance_specific_mutation_map <- sample_resistance_mutation[which(!is.na(match(sample_resistance_mutation[,3],specific_mutation))),]
	common_mutation <- merge(sample_resistance_mutation,sensitive_mutation,by="ID")
	if(dim(common_mutation)[1]>0){
		for(ID in unique(common_mutation[,1])){
			each_ID <- common_mutation[which(common_mutation[,1]==ID),]
			each_ID$diff_freq <- mean(as.numeric(each_ID[,4]))/mean(as.numeric(each_ID[,7]))
			common_diff_mutation <- each_ID[which(each_ID$diff_freq>2),c(2,3,4,1)]
			colnames(common_diff_mutation) <- colnames(sample_resistance_mutation)
			diff_mutation <- rbind(resistance_specific_mutation_map,common_diff_mutation)
		}
	}
	else{
		diff_mutation <- resistance_specific_mutation_map
	}
	return(diff_mutation)
}