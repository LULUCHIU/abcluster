#!/home/Tina/miniconda3/bin/R

#' @create by Jan 10 2022
#' @author by Tina


library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(MASS)
library(gdata)
library(argparser)

# argv <- arg_parser("")
# argv <- add_argument(argv,"--infile", help="")
# 
# argv <- parse_args(argv)
# 
# infile<-as.character(argv$infile)
AnnotateCDR<-function(myData,pyroute,dirpath){

infile<-myData



VH<-infile[,c("ID","VH")]
VL<-infile[,c("ID","VL")]

VHinput<-VH$VH
VLinput<-VL$VL



for(vh in 1:length(VHinput)){
  
if (is.na(VHinput[vh])){
   data_all_Hm<-data.frame(H_FW1=NA,H_CDR1=NA,H_FW2=NA,H_CDR2=NA,H_FW3=NA,H_CDR3=NA,H_FW4=NA,ID=VH$ID[vh])
   if(file.exists(paste0(dirpath,"/data_use_H.txt"))!="TRUE"){
  data_em <- data.frame(matrix(ncol = 8, nrow = 0))
  x <- c("H_FW1", "H_CDR1", "H_FW2","H_CDR2","H_FW3","H_CDR3","H_FW4","ID")
  colnames(data_em) <- x
  write.table(data_em,paste0(dirpath,"/data_use_H.txt"),quote = F,sep="\t")
  } 

  data_em<-read.table(paste0(dirpath,"/data_use_H.txt"),header=T,sep="\t")

}else{

system(paste0('/home/shiny/miniconda3/bin/python ',pyroute,'/run.py',' ','"',VHinput[vh],'"',' ',dirpath),TRUE)

  
if(file.exists(paste0(dirpath,"/data_use_H.txt"))!="TRUE"){
  data_em <- data.frame(matrix(ncol = 8, nrow = 0))
  x <- c("H_FW1", "H_CDR1", "H_FW2","H_CDR2","H_FW3","H_CDR3","H_FW4","ID")
  colnames(data_em) <- x
  write.table(data_em,paste0(dirpath,"/data_use_H.txt"),quote = F,sep="\t")
}

data_em<-read.table(paste0(dirpath,"/data_use_H.txt"),header=T,sep="\t")
t1<-NA
try({t1<-read.table(paste0(dirpath,"/test.txt"),header=F,fill=TRUE,sep=" ",na.strings ="", stringsAsFactors= F)},silent=TRUE)
#t1<-read.table(paste0(dirpath,"/test.txt"),header=F,sep=" ",fill=TRUE,na.strings ="", stringsAsFactors= F)
if(length(t1)==1 & is.na(t1)){
  data_all_Hm<-data.frame(H_FW1="WrongSeq",H_CDR1=NA,H_FW2=NA,H_CDR2=NA,H_FW3=NA,H_CDR3=NA,H_FW4=NA,ID=VH$ID[vh])
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


data_all_Hm<-data.frame(H_FW1=fw1,H_CDR1=cdr1,H_FW2=fw2,H_CDR2=cdr2,H_FW3=fw3,H_CDR3=cdr3,H_FW4=fw4,ID=VH$ID[vh])
}
}
data_all_H<-data.frame(rbind(data_em,data_all_Hm))
write.table(data_all_H,paste0(dirpath,"/data_use_H.txt"),sep="\t",quote=F,row.names=F)
#write.table(data_all,"data_use_H.txt",sep="\t",quote = F)

}
 

  for(vl in 1:length(VLinput)){
 	if (is.na(VLinput[vl])){
	data_all_Lm<-data.frame(L_FW1=NA,L_CDR1=NA,L_FW2=NA,L_CDR2=NA,L_FW3=NA,L_CDR3=NA,L_FW4=NA,ID=VL$ID[vl])
	if(file.exists(paste0(dirpath,"/data_use_L.txt"))!="TRUE"){
      data_em <- data.frame(matrix(ncol = 8, nrow = 0))
      x <- c("L_FW1", "L_CDR1", "L_FW2","L_CDR2","L_FW3","L_CDR3","L_FW4","ID")
      colnames(data_em) <- x
      write.table(data_em,paste0(dirpath,"/data_use_L.txt"),quote = F,sep="\t")
    }
    data_em<-read.table(paste0(dirpath,"/data_use_L.txt"),header=T,sep="\t")	

  }else{	   
    #system(paste0("python run.py"," ",VLinput[vl]),TRUE)
    #system(paste0("python ",pyroute,"/run.py"," ",VLinput[vl]),TRUE)
    system(paste0('/home/shiny/miniconda3/bin/python ',pyroute,'/run.py',' ','"',VLinput[vl],'"',' ',dirpath),TRUE)    

    
    if(file.exists(paste0(dirpath,"/data_use_L.txt"))!="TRUE"){
      data_em <- data.frame(matrix(ncol = 8, nrow = 0))
      x <- c("L_FW1", "L_CDR1", "L_FW2","L_CDR2","L_FW3","L_CDR3","L_FW4","ID")
      colnames(data_em) <- x
      write.table(data_em,paste0(dirpath,"/data_use_L.txt"),quote = F,sep="\t")
    }
    
    data_em<-read.table(paste0(dirpath,"/data_use_L.txt"),header=T,sep="\t")
    t1<-NA
    try({t1<-read.table(paste0(dirpath,"/test.txt"),header=F,fill=TRUE,sep=" ",na.strings ="", stringsAsFactors= F)},silent=TRUE)   
# t1<-read.table(paste0(dirpath,"/test.txt"),header=F,sep=" ",fill=TRUE,na.strings ="", stringsAsFactors= F)
   if(length(t1)==1 & is.na(t1)){
     data_all_Lm<-data.frame(L_FW1="WrongSeq",L_CDR1=NA,L_FW2=NA,L_CDR2=NA,L_FW3=NA,L_CDR3=NA,L_FW4=NA,ID=VL$ID[vl])
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
    try({fw1<-paste(c(as.matrix(t1[1,c(1:(cdr1posi[1]-1))])),collapse = "")},silent=TRUE)
   # if(!exists("fw1")){
   # fw1<-NA
   # }
    fw2<-NA
    try({fw2<-paste(c(as.matrix(t1[1,c((cdr1posi[length(cdr1posi)]+1):(cdr2posi[1]-1))])),collapse = "")},silent=TRUE)
   # if(!exists("fw2")){
   # fw2<-NA
   # }
    fw3<-NA
    try({fw3<-paste(c(as.matrix(t1[1,c((cdr2posi[length(cdr2posi)]+1):(cdr3posi[1]-1))])),collapse = "")},silent=TRUE)
   # if(!exists("fw3")){
   # fw3<-NA
   # }
    fw4<-NA
    try({fw4<-paste(c(as.matrix(t1[1,c((cdr3posi[length(cdr3posi)]+1):dim(t1)[2])])),collapse = "")},silent=TRUE)
   # if(!exists("fw4")){
   # fw4<-NA
   # }
    data_all_Lm<-data.frame(L_FW1=fw1,L_CDR1=cdr1,L_FW2=fw2,L_CDR2=cdr2,L_FW3=fw3,L_CDR3=cdr3,L_FW4=fw4,ID=VL$ID[vl])
  }
  }
   data_all_L<-data.frame(rbind(data_em,data_all_Lm))
    write.table(data_all_L,paste0(dirpath,"/data_use_L.txt"),sep="\t",quote = F,row.names=F)
  }
  #data_all_L<-cbind(VL,data_all_L)
  #write.table(data_all_L,"data_use_L.txt",sep="\t",quote = F)
  data_all_H<-read.table(paste0(dirpath,"/data_use_H.txt"),sep="\t",header=T)
  data_all_L<-read.table(paste0(dirpath,"/data_use_L.txt"),sep="\t",header=T)

data_all_H<-data_all_H[,c("ID","H_FW1", "H_CDR1", "H_FW2","H_CDR2","H_FW3","H_CDR3","H_FW4")]
data_all_L<-data_all_L[,c("ID","L_FW1", "L_CDR1", "L_FW2","L_CDR2","L_FW3","L_CDR3","L_FW4")]

  data_all<-merge(infile,data_all_H,by="ID")
  data_all<-merge(data_all,data_all_L,by="ID")
  write.table(data_all,paste0(dirpath,"/data_output.txt"),sep="\t",quote = F,row.names=F)
}
