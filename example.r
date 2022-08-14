patient_list <- read.table("patient_list.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
LGG_patient <- patient_list[which(patient_list[,1]=="LGG_Temozolomide"),]
resistant_mutation <- read.table("LGG_Temozolomide_resistant_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
sensitive_mutation <- read.table("LGG_Temozolomide_sensitive_mutation.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
resistant_exp <- read.table("LGG_Temozolomide_resistant_exp.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
sensitive_exp <- read.table("LGG_Temozolomide_sensitive_exp.txt",header=T,sep="\t",quote="",stringsAsFactors=F,row.names=1)
GRN <- read.table("GRN.txt",header=T,sep="\t",quote="",stringsAsFactors=F)
all_driver <- c()
for(sample_ID in unique(LGG_patient[,2])){
    each_result <- DRdriver(sample_ID)
	all_driver <- rbind(all_driver,each_result)
}
write.table(all_result,"driver.txt",col.names=T,row.names=F,sep="\t",quote=F)