#!/home/shiny/miniconda3/bin/R
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(argparser))

argv<-commandArgs(TRUE)
infile<-argv[1]
dirpath<-argv[2]
outfile<-argv[3]
pyroute<-argv[4]


#myData<-VH_tmp

#outfile<-"data_use_H.txt"

#AnnotateCDR<-function(myData,pyroute,dirpath,outfile){
  #infile<-infile
  infile<-read.table(infile,header=T,sep="\t")
  #infile<-read.table("VH_test.txt",header=T,sep="\t")
  #outfile<-"H_output.txt"
  #pyroute<-"/home/luting/script"
  #dirpath<-"/home/luting/EMBOSS_test"
  #inputseq<-infile[,2]
  
  for(vh in 1:length(infile[,2])){
    
    if (is.na(infile[,2][vh])){
      data_all_Hm<-data.frame(FW1=NA,CDR1=NA,FW2=NA,CDR2=NA,FW3=NA,CDR3=NA,FW4=NA,ID=infile$ID[vh])
      if(file.exists(paste0(dirpath,"/",outfile))!="TRUE"){
        data_em <- data.frame(matrix(ncol = 8, nrow = 0))
        x <- c("FW1", "CDR1", "FW2","CDR2","FW3","CDR3","FW4","ID")
        colnames(data_em) <- x
        write.table(data_em,paste0(dirpath,"/",outfile),quote = F,sep="\t")
      } 
      
      data_em<-read.table(paste0(dirpath,"/",outfile),header=T,sep="\t")
      
    }else{
      
      system(paste0('/home/shiny/miniconda3/bin/python ',pyroute,'/run.py',' ','"',infile[,2][vh],'"',' ',dirpath),TRUE)
      
      
      if(file.exists(paste0(dirpath,"/",outfile))!="TRUE"){
        data_em <- data.frame(matrix(ncol = 8, nrow = 0))
        x <- c("FW1", "CDR1", "FW2","CDR2","FW3","CDR3","FW4","ID")
        colnames(data_em) <- x
        write.table(data_em,paste0(dirpath,"/",outfile),quote = F,sep="\t")
      }
      
      data_em<-read.table(paste0(dirpath,"/",outfile),header=T,sep="\t")
      t1<-NA
      try({t1<-read.table(paste0(dirpath,"/test.txt"),header=F,fill=TRUE,sep=" ",na.strings ="", stringsAsFactors= F)},silent=TRUE)
      #t1<-read.table(paste0(dirpath,"/test.txt"),header=F,sep=" ",fill=TRUE,na.strings ="", stringsAsFactors= F)
      if(length(t1)==1 & is.na(t1)){
        data_all_Hm<-data.frame(FW1="WrongSeq",CDR1=NA,FW2=NA,CDR2=NA,FW3=NA,CDR3=NA,FW4=NA,ID=infile$ID[vh])
      }else{
        t1<-data.frame(t(t1))
        
        #H10-
        posi<-which(!is.na(t1[7,]))
        t1[c(8:9),posi]<-t1[c(7:8),posi]
        
        
        #CDR1,2
        posi2<-which(!is.na(t1[2,]))
        t1[c(7:9),posi2]<-t1[c(2:4),posi2]
        
        #H100-
        posi1<-which(!is.na(t1[6,]))
        t1[c(8:9),posi1]<-t1[c(6:7),posi1]
        
        #CDR3
        posi3<-which(!is.na(t1[1,]))
        t1[c(7:9),posi3]<-t1[c(1:3),posi3]
        
        
        t1<-t1[c(7:9),]
        colnames(t1)<-t1[1,]
        colnames(t1)<-rep("fr",length(colnames(t1)))
        
        colnames(t1)[which(t1[1,]=="CDR1")]<-"CDR1"
        colnames(t1)[which(t1[1,]=="CDR2")]<-"CDR2"
        colnames(t1)[which(t1[1,]=="CDR3")]<-"CDR3"
        #t1[1,which(!t1[1,] %in% c("CDR1","CDR2","CDR3"))]<-"fr"
        #t1[1,c(which(!grepl("CDR",t1[1,])))]<-rep("fr",length(c(which(!grepl("CDR",t1[1,])))))
        #colnames(t1)<-t1[1,]
        t1<-t1[-c(1,2),]
        
        
        cdr1posi<-which(colnames(t1)=="CDR1")
        if(length(cdr1posi)!=0){
          cdr1<-paste(c(as.matrix(t1[1,cdr1posi])),collapse = "")
        }else{
          cdr1<-NA
        }
        cdr2posi<-which(colnames(t1)=="CDR2")
        if(length(cdr2posi)!=0){
          cdr2<-paste(c(as.matrix(t1[1,cdr2posi])),collapse = "")
        }else{
          cdr2<-NA
        }
        cdr3posi<-which(colnames(t1)=="CDR3")
        if(length(cdr3posi)!=0){
          cdr3<-paste(c(as.matrix(t1[1,cdr3posi])),collapse = "")
        }else{
          cdr3<-NA
        }
        
        fw1<-NA
        try({fw1<-paste(c(as.matrix(t1[1,c(1:cdr1posi[1]-1)])),collapse = "")},silent=TRUE)
        #if(!exists("fw1")){
        #fw1<-NA
        #}
        
        fw2<-NA
        try({fw2<-paste(c(as.matrix(t1[1,c((cdr1posi[length(cdr1posi)]+1):(cdr2posi[1]-1))])),collapse = "")},silent=TRUE)
        #if(!exists("fw2")){
        #fw2<-NA
        #}
        
        fw3<-NA
        try({fw3<-paste(c(as.matrix(t1[1,c((cdr2posi[length(cdr2posi)]+1):(cdr3posi[1]-1))])),collapse = "")},silent=TRUE)
        #if(!exists("fw3")){
        #fw3<-NA
        #}
        
        fw4<-NA
        try({fw4<-paste(c(as.matrix(t1[1,c((cdr3posi[length(cdr3posi)]+1):dim(t1)[2])])),collapse = "")},silent=TRUE)
        #if(!exists("fw4")){
        #fw4<-NA
        #}
        
        
        data_all_Hm<-data.frame(FW1=fw1,CDR1=cdr1,FW2=fw2,CDR2=cdr2,FW3=fw3,CDR3=cdr3,FW4=fw4,ID=infile$ID[vh])
      }
    }
    data_all_H<-data.frame(rbind(data_em,data_all_Hm))
    write.table(data_all_H,paste0(dirpath,"/",outfile),sep="\t",quote=F,row.names=F)
    #write.table(data_all,"data_use_H.txt",sep="\t",quote = F)
    
  }
  
#}
