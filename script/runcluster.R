#!/home/Tina/miniconda3/bin/R
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(MASS)
library(gdata)

 #Heavy Chain Framework seq (A0A0C4DH42 (HV366_HUMAN))
  hfr1<-"EVQLVESGGGLIQPGGSLRLSCAASGFTVS"
  hfr2<-"WVRQAPGKGLEWVS"
  hfr3<-"RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
  #Light Chain Framework seq (>pdb|3FZU|L Chain L, immunoglobulin IgG1 Fab, light chain)
  lfr1<-"DIQMTQSPSSLSASVGDRVTITC"
  lfr2<-"WYQQKPGKAPKLLIY"
  lfr3<-"GVPSRFSGSGSGTEFTLTVSSLQPEDFATYYC"

  #infile<-read.xlsx("/Users/qiuluting/Desktop/CTMbio/Needs/Shinyapp/abcluster/Testfile/yourname@Mode1_example.xlsx")
  #infile<-read.xlsx("/data/shiny/abcluster/luting_qiu/Mode1_example_Tue_Mar__1_10:20:16_CST_2022/luting_qiu@Mode1_example.xlsx")
  #outdir<-"/data/shiny/abcluster/luting_qiu/Mode1_example_Tue_Mar__1_10:20:16_CST_2022"
#==============HC only
runcluster<-function(infile,outdir,seqid,cover,mode){
if(mode %in% c("vh","vhvlcdr")){

data_ori<-infile
data<-data_ori[,c("ID","H_CDR1","H_CDR2","H_CDR3")]
data<-data[order(data$H_CDR3,data$H_CDR1,data$H_CDR2),]
#sort CDR3 CDR1 CDR2
data$Full<-paste0(data$H_CDR1,data$H_CDR2,data$H_CDR3)
#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
##=========================================
if(dim(dupset)[1]!=0){
countdup<-data.frame(table(dupset$Full))
names(countdup)[1]<-"Full"
dupset<-merge(dupset,countdup,by="Full",all.x=T)
countdup<-dupset[,c("ID","Freq")]
}
##==========================================
if(length(which(duplicated(data$Full))!=0)){
  overlapid<-data$ID[which(duplicated(data$Full))]
  data1<-data[-which(duplicated(data$Full)),]
  retainid<-intersect(dupset$ID,data1$ID)
  countdup$dupID<-rep(NA,dim(countdup)[1])
  for(dd in 1:length(retainid)){
  dupid<-setdiff(dupset$ID[which(dupset$Full==dupset$Full[which(dupset$ID==retainid[dd])])],retainid[dd])
  countdup$dupID[which(countdup$ID %in% dupid)]<-retainid[dd]  
  }
}else{
  data1<-data
  overlapid<-c() 
}

##=========================================


data$Name<-paste0(">",data$ID)

data$Com<-paste0(data$Name,"\n",data$H_CDR3,hfr3,data$H_CDR2,hfr2,data$H_CDR1)
out<-data.frame(data$Com)
write.table(out,paste0(outdir,"/testHC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/testHC.fasta"," ",outdir,"/result_HC ",outdir,"/tmp ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_HC_cluster.tsv > ',outdir,'/result_HC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_HC_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")
parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"H_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")
output<-merge(data_ori,data.frame(output[,c("Member","H_cluster")]),by.x="ID",by.y="Member")
output$H_cluster<-reorder.factor(output$H_cluster,new.order=datadb$H_cluster)
output<-output %>% arrange(H_cluster)

if(length(overlapid)!=0){
  dupset<-output[which(output$ID %in% dupset$ID),]
  dupset<-merge(dupset,countdup,by="ID",all.x=T)
  dupset$Retain<-ifelse(dupset$ID %in% retainid,"O","X")
  write.table(dupset,paste0(outdir,"/duplicated_seqs_HC.txt"),sep="\t",col.names = T,row.names = F,quote=F)
  output<-output[-which(output$ID %in% overlapid),]
  write.xlsx(output,paste0(outdir,"/Cluster_HC_Result.xlsx"))
}else{
   write.xlsx(output,paste0(outdir,"/Cluster_HC_Result.xlsx"))
  }
  #write.xlsx(output,paste0(outdir,"/Cluster_HC_Result.xlsx"))
  if (mode=="vh"){
  Outlist<-list()
  #names(Outlist)<-c("Cluster_HC_Result","Duplicated_set")
  Outlist[[1]]<-output
  if(length(overlapid)!=0){
  Outlist[[2]]<-dupset
 }else{
  Outlist[[2]]<-data.frame(Noduplicated=c("Noduplicated"))
 }
  names(Outlist)<-c("Cluster_HC_Result","Duplicated_set")
  write.xlsx(Outlist,paste0(outdir,"/Cluster_Output.xlsx"))
                 }
#HC_LC===============================
if(mode=="vhvlcdr"){
data_ori<-infile

data<-data_ori[,c("ID","H_CDR1","H_CDR2","H_CDR3","L_CDR1","L_CDR2","L_CDR3")]
data<-data[order(data$H_CDR3,data$H_CDR1,data$H_CDR2,data$L_CDR3,data$L_CDR1,data$L_CDR2),]
#sort CDR3 CDR1 CDR2
data$Full<-paste0(data$H_CDR1,data$H_CDR2,data$H_CDR3,data$L_CDR1,data$L_CDR2,data$L_CDR3)

#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
##=========================================
if(dim(dupset)[1]!=0){
  countdup<-data.frame(table(dupset$Full))
  names(countdup)[1]<-"Full"
  dupset<-merge(dupset,countdup,by="Full",all.x=T)
  countdup<-dupset[,c("ID","Freq")]
}
##==========================================
if(length(which(duplicated(data$Full))!=0)){
  overlapid<-data$ID[which(duplicated(data$Full))]
  data1<-data[-which(duplicated(data$Full)),]
  retainid<-intersect(dupset$ID,data1$ID)
  countdup$dupID<-rep(NA,dim(countdup)[1])
  for(dd in 1:length(retainid)){
    dupid<-setdiff(dupset$ID[which(dupset$Full==dupset$Full[which(dupset$ID==retainid[dd])])],retainid[dd])
    countdup$dupID[which(countdup$ID %in% dupid)]<-retainid[dd]  
  }
}else{
  overlapid<-c()
  data1<-data
}

data$Name<-paste0(">",data$ID)
data$Com<-paste0(data$Name,"\n",data$H_CDR3,hfr3,data$H_CDR2,hfr2,data$H_CDR1,hfr1,data$L_CDR3,lfr3,
                 data$L_CDR2,lfr2,data$L_CDR1,lfr1)
out<-data.frame(data$Com)

write.table(out,paste0(outdir,"/test_HCLC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/test_HCLC.fasta"," ",outdir,"/result_HCLC ",outdir,"/tmp_HCLC ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_HCLC_cluster.tsv > ',outdir,'/result_HCLC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_HCLC_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")
parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"HC_LC_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")
output<-merge(data_ori,data.frame(output[,c("Member","HC_LC_cluster")]),by.x="ID",by.y="Member")
output$HC_LC_cluster<-reorder.factor(output$HC_LC_cluster,new.order=datadb$HC_LC_cluster)
output<-output %>% arrange(HC_LC_cluster)

if(length(overlapid)!=0){
  dupset<-output[which(output$ID %in% dupset$ID),]
  #
  dupset<-merge(dupset,countdup,by="ID",all.x=T)
  dupset$Retain<-ifelse(dupset$ID %in% retainid,"O","X")
  write.table(dupset,paste0(outdir,"/duplicated_seqs_HCLC.txt"),sep="\t",col.names = T,row.names = F,quote=F)
  output<-output[-which(output$ID %in% overlapid),]
}

write.xlsx(output,paste0(outdir,"/Cluster_HCLC_Result.xlsx"))

##LC===========================

data_ori<-infile

data<-data_ori[,c("ID","L_CDR1","L_CDR2","L_CDR3")]
data<-data[order(data$L_CDR3,data$L_CDR1,data$L_CDR2),]
#sort CDR3 CDR1 CDR2
data$Full<-paste0(data$L_CDR1,data$L_CDR2,data$L_CDR3)

#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]

##=========================================
if(dim(dupset)[1]!=0){
  countdup<-data.frame(table(dupset$Full))
  names(countdup)[1]<-"Full"
  dupset<-merge(dupset,countdup,by="Full",all.x=T)
  countdup<-dupset[,c("ID","Freq")]
}

##==========================================
if(length(which(duplicated(data$Full))!=0)){
  overlapid<-data$ID[which(duplicated(data$Full))]
  data1<-data[-which(duplicated(data$Full)),]
  retainid<-intersect(dupset$ID,data1$ID)
  countdup$dupID<-rep(NA,dim(countdup)[1])
  for(dd in 1:length(retainid)){
    dupid<-setdiff(dupset$ID[which(dupset$Full==dupset$Full[which(dupset$ID==retainid[dd])])],retainid[dd])
    countdup$dupID[which(countdup$ID %in% dupid)]<-retainid[dd]  
  }
}else{
  overlapid<-c()
  data1<-data
}

##=========================================
data$Name<-paste0(">",data$ID)

data$Com<-paste0(data$Name,"\n",data$L_CDR3,lfr3,data$L_CDR2,lfr2,data$L_CDR1)
out<-data.frame(data$Com)
write.table(out,paste0(outdir,"/testLC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
#write.table(dupset,paste0(outdir,"/duplicated_seqs_LC.txt"),sep="\t",col.names = F,row.names = F,quote=F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/testLC.fasta"," ",outdir,"/result_LC ",outdir,"/tmp_LC ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_LC_cluster.tsv > ',outdir,'/result_LC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_LC_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")
parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"L_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")
output<-merge(data_ori,data.frame(output[,c("Member","L_cluster")]),by.x="ID",by.y="Member")
output$L_cluster<-reorder.factor(output$L_cluster,new.order=datadb$L_cluster)
output<-output %>% arrange(L_cluster)

if(length(overlapid)!=0){
  dupset<-output[which(output$ID %in% dupset$ID),]
  #
  dupset<-merge(dupset,countdup,by="ID",all.x=T)
  dupset$Retain<-ifelse(dupset$ID %in% retainid,"O","X")
  write.table(dupset,paste0(outdir,"/duplicated_seqs_LC.txt"),sep="\t",col.names = T,row.names = F,quote=F)
  output<-output[-which(output$ID %in% overlapid),]
}
write.xlsx(output,paste0(outdir,"/Cluster_LC_Result.xlsx"))

}
}
}


