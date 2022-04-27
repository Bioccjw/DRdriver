##@@=====================================================================================================
##@@DRdriver: identifying drug resistance driver genes using individual-specific gene regulatory network
##@@genetic algorithm for identifying drug resistance driver genes

##@@network: the individual-specific network with three columns
##first column: candidate driver genes
##second column: corresponding targets
##third column: whether the targets were DEGs. Only "diff" and "non-diff" were promised.

##@@score: all candidates
##@@NUMPOP: the number of individuals in the population
##@@STARTNODENUM: the proportion of genes that retained from old population
##@@SELECTRATE: the proportion of selection
##@@CROSSOVERRATE: crossover rate
##@@VARIATIONRATE: mutation rate
##@@ITERATION: the time of iteration

Genetic_Algorithm <- function(network,
   score,
   NUMPOP = 1000,
   STARTNODENUM = 0.5,
   SELECTRATE = 0.5,
   CROSSOVERRATE = 0.9,
   VARIATIONRATE = 0.3,
   ITERATION = 100){
		
   ##compute the fitness
   Fitness <- function(pop,network,score){ 
        scoreCompute <- function(X){
             sel_position <- which(X==1)
             sel_node <- score[sel_position,1]
             sel_control <- network[which(!is.na(match(network[,1],sel_node))),]
             desire_score <- length(unique(sel_control[which(sel_control[,3]=="diff"),2]))/length(unique(network[which(network[,3]=="diff"),2]))
             undesire_score <- length(unique(sel_control[which(sel_control[,3]=="non-diff"),2]))/length(unique(network[which(network[,3]=="non-diff"),2]))
             final_score <- desire_score-undesire_score
             return(final_score)
        }	
        fitness <- apply(pop,1,scoreCompute)
        return(fitness)
   }
   ##select Elite
   Elite <- function(pop,fitness,SELECTRATE){
        order_number <- order(fitness,decreasing=T)
        sel_parent_number <- ceiling(nrow(pop)*SELECTRATE)
        Elite <- pop[order_number[1:sel_parent_number],]
        return(Elite)
   }

   ParentPop <- function(pop,fitness){
        normalized_fitness <- (fitness-min(fitness))/(max(fitness)-min(fitness))
        sumFitness <- sum(normalized_fitness)
        accP_1 <- normalized_fitness/sumFitness
        accP_1 <- accP_1[order(accP_1,decreasing=T)]
        pop_1 <- pop[order(accP_1,decreasing=T),]
        accP <- cumsum(accP_1)
        sel_parent_number <- ceiling(nrow(pop)*SELECTRATE)
        ParentPop <- c()
        for(n in 1:sel_parent_number){
            each_sel <- pop_1[which(accP>=runif(1,min=0,max=1))[1],]
            ParentPop <- rbind(ParentPop,each_sel)
        }
        return(ParentPop)
   }

   ##crossover
   Crossover <- function(ParentPop,CROSSOVERRATE){
        X <- ParentPop[sample(nrow(ParentPop),1),]
        Y <- ParentPop[sample(nrow(ParentPop),1),]
        kid_1 <- X
        kid_2 <- Y
        if(runif(1) < CROSSOVERRATE){
            crossLocation <- ceiling((length(X)-1)*runif(1))+1
            kid_1[crossLocation:length(X)] <- Y[crossLocation:length(Y)]
            kid_2[crossLocation:length(Y)] <- X[crossLocation:length(X)]
        }
        kid_1 <- matrix(kid_1,ncol=LENGTH)
        kid_2 <- matrix(kid_2,ncol=LENGTH)
        kidsPop <- rbind(kid_1,kid_2)
        return(kidsPop)
   }

   ##Variation
   Variation <- function(ParentPop,VARIATIONRATE){
        X <- ParentPop[sample(nrow(ParentPop),1),]
        Y <- ParentPop[sample(nrow(ParentPop),1),]
        kid_1 <- X
        kid_2 <- Y
        if(runif(1) < VARIATIONRATE){
            location_1 <- ceiling((length(X)-1)*runif(1))+1
            location_2 <- ceiling((length(X)-1)*runif(1))+1
            kid_1[location_1] <- 1-kid_1[location_1]
            kid_2[location_2] <- 1-kid_2[location_2]
        }
        kid_1 <- matrix(kid_1,ncol=LENGTH)
        kid_2 <- matrix(kid_2,ncol=LENGTH)
        kidsPop <- rbind(kid_1,kid_2)
        return(kidsPop)
   }

   ##Incoding
   Incoding <- function(best_pop,score){
        sel_node <- which(best_pop==1)
        node <- score[sel_node,1]
        return(node)
   }

   ##deleting
   delete_gene <- function(X){
        sel_position <- sample(length(which(X==1)),ceiling(length(which(X==1))*(1-STARTNODENUM_1)))
        X[which(X==1)[sel_position]] <- 0
        return(X)
   }



   LENGTH <- length(score[,1])
   first_pop <- matrix(rep(1,LENGTH*NUMPOP),ncol=LENGTH)
   first_fitness <- Fitness(first_pop,network,score)
   sel_best_pop <- first_pop
   sel_best_fitness <- first_fitness
   best_fitness <- 0
   best_result <- score[,1]
   all_fitness <- c()
   gene_number <- c()
   zero_time <- 0
   STARTNODENUM_1 <- STARTNODENUM
   for(time in 1:ITERATION){
        if(zero_time>=5){
            STARTNODENUM_1 <- 1
        }
        pop <- t(apply(sel_best_pop,1,delete_gene))
        fitness <- Fitness(pop,network,score)
        elite <- Elite(pop,fitness,SELECTRATE)
        parentpop <- ParentPop(pop,fitness)
        kidsPop <- c()
        n=1
        while(n <= (nrow(pop)-nrow(parentpop))){
            each_temp_pop <- c()
            for(i in 1:10){
                kidsPop_1 <- Crossover(parentpop,CROSSOVERRATE)
                kidsPop_2 <- Variation(kidsPop_1,VARIATIONRATE)
                each_temp_pop <- rbind(each_temp_pop,kidsPop_2)
            }
            each_temp_fitness <- Fitness(each_temp_pop,network,score)
            sel_fitness_position <- which(each_temp_fitness==max(each_temp_fitness))
            sel_kidpop <- each_temp_pop[sel_fitness_position[1],]
            kidsPop <- rbind(kidsPop,sel_kidpop)
            n <- n+1
            #print(n)
        }
        temp_pop <- rbind(elite,kidsPop)
        new_fitness <- Fitness(temp_pop,network,score)
        delta_fitness <- max(new_fitness)-best_fitness
        print(delta_fitness)
        if(delta_fitness>0){
            sel_best_pop <- temp_pop
            sel_best_fitness <- new_fitness
        }
        best_position <- order(sel_best_fitness,decreasing=T)[1]
        best_pop_1 <- sel_best_pop[best_position,]
        best_fitness <- sel_best_fitness[best_position]
        best_pop <- matrix(best_pop_1,ncol=LENGTH)
        best_result <- Incoding(best_pop,score)
        print(length(best_result))
        all_fitness <- c(all_fitness,best_fitness)
        gene_number <- c(gene_number,length(best_result))
        pop <- sel_best_pop		
        if(delta_fitness<=0){
            zero_time <- zero_time+1
        }
        else{
            zero_time <- 0
        }
        if(zero_time >= 10){
            break
        }
        print(zero_time)
        print(time)	
        }
		
        result <- list(best_result,best_fitness,time,all_fitness,gene_number)
        return(result)
}
