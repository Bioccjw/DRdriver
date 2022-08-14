##the function was used for transform the sample ID
##the sample ID contains "-". When reading them as column names in exp.txt, the "-" always become ".".
sample_transformation <- function(x){
    split_name <- unlist(strsplit(x,split="-"))
	merge_1 <- paste(split_name[1],split_name[2],sep=".")
	merge_2 <- paste(merge_1,split_name[3],sep=".")
	return(merge_2)
}
patient_list <- read.table("patient_list.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
patient_list[,4] <- unlist(lapply(patient_list[,2],sample_transformation))
LGG_patient <- patient_list[which(patient_list[,1]=="LGG_Temozolomide"),]
resistant_mutation <- read.table("LGG_Temozolomide_resistant_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
sensitive_mutation <- read.table("LGG_Temozolomide_sensitive_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
resistant_EXP <- read.table("LGG_Temozolomide_resistant_exp.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
sensitive_EXP <- read.table("LGG_Temozolomide_sensitive_exp.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
GRN <- read.table("GRN.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
all_driver <- c()
for(sample_ID in unique(resistant_patient[,2])){
    each_result <- DRdriver(sample_ID)
	all_driver <- rbind(all_driver,each_result)
}
write.table(all_result,"driver.txt",col.names=T,row.names=F,sep="\t",quote=F)