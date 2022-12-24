#!/home/Tina/miniconda3/bin/R
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(MASS)
library(gdata)

 #Heavy Chain Framework seq (A0A0C4DH42 (HV366_HUMAN))
  #hfr1<-"EVQLVESGGGLIQPGGSLRLSCAASGFTVS"
  #hfr2<-"WVRQAPGKGLEWVS"
  #hfr3<-"RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
 #Light Chain Framework seq (>pdb|3FZU|L Chain L, immunoglobulin IgG1 Fab, light chain)
  #lfr1<-"DIQMTQSPSSLSASVGDRVTITC"
  #lfr2<-"WYQQKPGKAPKLLIY"
  #lfr3<-"GVPSRFSGSGSGTEFTLTVSSLQPEDFATYYC"
 #ori_fw
  fr1<-"EIVLTQSPGTLSLSPGERATLSC"
  fr2<-"IGWFRQAPGKEREGVAW"
  fr3<-"YYADSVKGRFTISSDNAKNTVSLQMNSLKPEDTAKYYC" 

#==============HC only
runcluster<-function(infile,outdir,seqid,cover){

data_ori<-infile
data<-data_ori[,c("ID","H_CDR1","H_CDR2","H_CDR3")]
data<-data[order(data$H_CDR3,data$H_CDR1,data$H_CDR2),]
#sort CDR3 CDR1 CDR2
data$Full<-paste0(data$H_CDR1,data$H_CDR2,data$H_CDR3)

#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
#data<-data[-which(duplicated(data$Full)),]
if(length(which(duplicated(data$Full))!=0)){
  data<-data[-which(duplicated(data$Full)),]
}else{
  data<-data
}
data$Name<-paste0(">",data$ID)
#data$Com<-paste0(data$Name,"\n",data$CDR1,fr2,data$CDR2,fr3,data$CDR3)
data$Com<-paste0(data$Name,"\n",data$H_CDR3,fr3,data$H_CDR2,fr2,data$H_CDR1)
out<-data.frame(data$Com)
write.table(out,paste0(outdir,"/testHC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
write.table(dupset,paste0(outdir,"/duplicated_seqs_HC.txt"),sep="\t",col.names = F,row.names = F,quote=F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/testHC.fasta"," ",outdir,"/result_HC ",outdir,"/tmp ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_HC_cluster.tsv > ',outdir,'/result_HC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_HC_cluster.txt"),sep="\t",header=F)
#output<-read.table(paste0(wdir,"/output_linux/result_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")
parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"H_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")
output$Represent<-gsub("_","-",output$Represent)
output$Member<-gsub("_","-",output$Member)
output<-merge(data_ori,data.frame(output[,c("Member","H_cluster")]),by.x="ID",by.y="Member")
output$H_cluster<-reorder.factor(output$H_cluster,new.order=datadb$H_cluster)
output<-output %>% arrange(H_cluster)

#write.table(output,paste0(outdir,"/Cluster_HC_Result.txt"),sep="\t",quote=F,row.names=F)
write.xlsx(output,paste0(outdir,"/Cluster_HC_Result.xlsx"))

#HC_LC===============================

#if (stringr::str_ends(infile,'txt')) {
#    data_ori <- read.table(infile, header = TRUE)
#    }else if (stringr::str_ends(infile,"(xlsx|xls)")) {
#      data_ori<-read.xlsx(infile)
#    }

data_ori<-infile

#data_ori<-read.xlsx(infile,sheet=1)
#data<-data_ori %>% select(ID,CDR1,CDR2,CDR3)
data<-data_ori[,c("ID","H_CDR1","H_CDR2","H_CDR3","L_CDR1","L_CDR2","L_CDR3")]
data<-data[order(data$H_CDR3,data$H_CDR1,data$H_CDR2,data$L_CDR3,data$L_CDR1,data$L_CDR2),]
#sort CDR3 CDR1 CDR2

data$Full<-paste0(data$H_CDR1,data$H_CDR2,data$H_CDR3,data$L_CDR1,data$L_CDR2,data$L_CDR3)

#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
if(length(which(duplicated(data$Full))!=0)){
  data<-data[-which(duplicated(data$Full)),]
}else{
  data<-data
}
data$Name<-paste0(">",data$ID)

data$Com<-paste0(data$Name,"\n",data$H_CDR3,fr3,data$H_CDR2,fr2,data$H_CDR1,fr1,data$L_CDR3,fr3,
                 data$L_CDR2,fr2,data$L_CDR1)
out<-data.frame(data$Com)


#write.table(out,paste0(wdir,"/test.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
write.table(out,paste0(outdir,"/test_HCLC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
write.table(dupset,paste0(outdir,"/duplicated_seqs_HCLC.txt"),sep="\t",col.names = F,row.names = F,quote=F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/test_HCLC.fasta"," ",outdir,"/result_HCLC ",outdir,"/tmp_HCLC ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_HCLC_cluster.tsv > ',outdir,'/result_HCLC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_HCLC_cluster.txt"),sep="\t",header=F)
#output<-read.table(paste0(wdir,"/output_linux/result_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")

parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"HC_LC_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")

output$Represent<-gsub("_","-",output$Represent)
output$Member<-gsub("_","-",output$Member)

output<-merge(data_ori,data.frame(output[,c("Member","HC_LC_cluster")]),by.x="ID",by.y="Member")

output$HC_LC_cluster<-reorder.factor(output$HC_LC_cluster,new.order=datadb$HC_LC_cluster)
output<-output %>% arrange(HC_LC_cluster)

#write.table(output,paste0(outdir,"/Cluster_HCLC_Result.txt"),sep="\t",quote=F,row.names=F)
write.xlsx(output,paste0(outdir,"/Cluster_HCLC_Result.xlsx"))

##LC===========================

#if (stringr::str_ends(infile,'txt')) {
#    data_ori <- read.table(infile, header = TRUE)
#    }else if (stringr::str_ends(infile,"(xlsx|xls)")) {
#      data_ori<-read.xlsx(infile)
#    }

data_ori<-infile
#data_ori<-read.xlsx(infile)
data<-data_ori[,c("ID","L_CDR1","L_CDR2","L_CDR3")]
data<-data[order(data$L_CDR3,data$L_CDR1,data$L_CDR2),]
#sort CDR3 CDR1 CDR2
data$Full<-paste0(data$L_CDR1,data$L_CDR2,data$L_CDR3)

#duplicated Full
dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
#data<-data[-which(duplicated(data$Full)),]
if(length(which(duplicated(data$Full))!=0)){
  data<-data[-which(duplicated(data$Full)),]
}else{
  data<-data
}
data$Name<-paste0(">",data$ID)
#data$Com<-paste0(data$Name,"\n",data$CDR1,fr2,data$CDR2,fr3,data$CDR3)
data$Com<-paste0(data$Name,"\n",data$L_CDR3,fr3,data$L_CDR2,fr2,data$L_CDR1)
out<-data.frame(data$Com)
write.table(out,paste0(outdir,"/testLC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
write.table(dupset,paste0(outdir,"/duplicated_seqs_LC.txt"),sep="\t",col.names = F,row.names = F,quote=F)

system(paste0("cd ",outdir),TRUE)
system(paste0("mmseqs easy-cluster ",outdir,"/testLC.fasta"," ",outdir,"/result_LC ",outdir,"/tmp_LC ","--min-seq-id ",seqid," -c ",
              cover," --cov-mode 0"," --cluster-reassign"," --alignment-mode 3"), TRUE)
system(paste0('sed -e "s/\r//g" ',outdir,'/result_LC_cluster.tsv > ',outdir,'/result_LC_cluster.txt'),TRUE)

#read output
output<-read.table(paste0(outdir,"/result_LC_cluster.txt"),sep="\t",header=F)
#output<-read.table(paste0(wdir,"/output_linux/result_cluster.txt"),sep="\t",header=F)
colnames(output)<-c("Represent","Member")
parent<-as.character(unique(output$Represent))
datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
colnames(datadb)[2]<-"L_cluster"
output<-merge(output,datadb,by.x="Represent",by.y="parent")
output$Represent<-gsub("_","-",output$Represent)
output$Member<-gsub("_","-",output$Member)
output<-merge(data_ori,data.frame(output[,c("Member","L_cluster")]),by.x="ID",by.y="Member")
output$L_cluster<-reorder.factor(output$L_cluster,new.order=datadb$L_cluster)
output<-output %>% arrange(L_cluster)

#write.table(output,paste0(outdir,"/Cluster_LC_Result.txt"),sep="\t",quote=F,row.names=F)
write.xlsx(output,paste0(outdir,"/Cluster_LC_Result.xlsx"))

}
#else{
#   #==============HC only
# if (stringr::str_ends(infile,'txt')) {
#     data_ori <- read.table(infile, header = TRUE)
#     }else if (stringr::str_ends(infile,"(xlsx|xls)")) {
#       data_ori<-read.xlsx(infile)
#     }
# 
#  
#  #data_ori<-read.xlsx(infile)
#   data<-data_ori[,c("ID","H_CDR1","H_CDR2","H_CDR3")]
#   data<-data[order(data$H_CDR3,data$H_CDR1,data$H_CDR2),]
#   #sort CDR3 CDR1 CDR2
#   data$Full<-paste0(data$H_CDR1,data$H_CDR2,data$H_CDR3)
#   
#   #duplicated Full
#   dupset<-data[data$Full %in% c(data$Full[which(duplicated(data$Full))]),]
#   #data<-data[-which(duplicated(data$Full)),]
#   if(length(which(duplicated(data$Full))!=0)){
#     data<-data[-which(duplicated(data$Full)),]
#   }else{
#     data<-data
#   }
#   data$Name<-paste0(">",data$ID)
#   #data$Com<-paste0(data$Name,"\n",data$CDR1,fr2,data$CDR2,fr3,data$CDR3)
#   data$Com<-paste0(data$Name,"\n",data$H_CDR3,fr3,data$H_CDR2,fr2,data$H_CDR1)
#   out<-data.frame(data$Com)
#   write.table(out,paste0(outdir,"/testHC.fasta"),sep="\t",col.names=F, row.names = F,quote = F)
#   write.table(dupset,paste0(outdir,"/duplicated_seqs_HC.txt"),sep="\t",col.names = F,row.names = F,quote=F)
#   
#   system(paste0("cd ",outdir),TRUE)
#   system(paste0("mmseqs easy-cluster ","testHC.fasta"," result_HC tmp ","--min-seq-id ",seqid," -c ",
#                 cover," --cov-mode ",mode," --cluster-reassign"," --alignment-mode 3"), TRUE)
#   system(paste0('sed -e "s/\r//g" ','result_HC_cluster.tsv > result_HC_cluster.txt'),TRUE)
#   
#   #read output
#   output<-read.table(paste0(outdir,"/result_HC_cluster.txt"),sep="\t",header=F)
#   #output<-read.table(paste0(wdir,"/output_linux/result_cluster.txt"),sep="\t",header=F)
#   colnames(output)<-c("Represent","Member")
#   parent<-as.character(unique(output$Represent))
#   datadb<-data.frame(cbind(parent,paste0("C",c(1:length(parent)))))
#   colnames(datadb)[2]<-"cluster"
#   output<-merge(output,datadb,by.x="Represent",by.y="parent")
#   output$Represent<-gsub("_","-",output$Represent)
#   output$Member<-gsub("_","-",output$Member)
#   output<-merge(data_ori,data.frame(output[,c("Member","cluster")]),by.x="ID",by.y="Member")
#   output$cluster<-reorder.factor(output$cluster,new.order=datadb$cluster)
#   output<-output %>% arrange(cluster)
#   
#   ######plot
#   # data<-output
#   # d1 <- data.frame(from="origin", to=paste("C", seq(1,length(unique(data$cluster))), sep=""))
#   # d2 <- data.frame(from=data$cluster,to=data$ID)
#   # edges <- rbind(d1, d2)
#   # 
#   # 
#   # # We can add a second data frame with information for each node!
#   # name <- unique(c(as.character(edges$from), as.character(edges$to)))
#   # vertices <- data.frame(
#   #   name=name,
#   #   group=c( rep(NA,(length(unique(data$cluster))+1)) , data$cluster)
#   # )
#   # # Create a graph object
#   # mygraph <- graph_from_data_frame(edges, vertices=vertices)
#   # 
#   # ggraph(mygraph, layout = 'dendrogram') + 
#   #   geom_edge_diagonal() +
#   #   geom_node_text(aes( label=name, filter=leaf, color=group) , angle=90 , hjust=1, nudge_y=-0.05,parse=TRUE,size=1.3) +
#   #   ylim(-.6, NA) +
#   #   theme(legend.position="none")
#   # ggsave(file=paste0(dirlist[ii],"_C",length(unique(data$cluster)),".pdf"),limitsize = FALSE,width=pretty(round(length(unique(data$cluster))))[2]/2.5,height =6)
#   write.table(output,paste0(outdir,"/Cluster_HC_Result.txt"),sep="\t",quote=F,row.names=F)
# }







