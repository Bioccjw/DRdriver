DRdriver <- function(sample_ID){
    network <- specific_network(sample_ID)
	score <- data.frame(unique(network[,1]))
	if(nrow(score)>1&length(unique(network[,3]))>1){
	    driver <- Genetic_Algorithm(network,score)
		return(driver)
	}
}