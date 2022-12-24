#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(DT)
library(ggraph)
library(igraph)
library(gdata)
library(dplyr)
library(openxlsx)
library(forcats)## for reordering the factor
library(scales)
library(gtools)
library(stringr)

#Sys.setenv(R_ZIPCMD="/usr/bin/zip")
ui<- shinyUI(fluidPage(
  theme=shinytheme("flatly"),
  navbarPage(
    "CtMbio-Tools",
    id = "main_navbar",
    
   # tabsetPanel(
    tabPanel(
      "Ab-Clustering",
      fluidRow(column(5,h3("Upload File")),
                column(6,h3("Result"))),
    #),
  #  ),
    
    sidebarLayout(
      sidebarPanel(
        helpText("This clustering method is based on MMseqs2, if you have interest in it
                 please browse their website by yourself : )"),
        #width=5,
        fileInput('file1', 'Choose Your File',
                  accept=c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv','.xlsx')),
        
        radioButtons("checkGroup","Which type you feed me?",
                           c("Mode1:VH+VL(CDR confirmed)"="vhvlcdr",
                             "Mode2:VH+VL(CDR NOT confirmed)"="vhvl",                                     
                             "Mode3:VH(CDR confirmed)"="vh")),
      
        sliderInput('condition_id',h4("Identity"),
                  min=0.5,max=1,step=0.01,value = 0.7),
       # sliderInput("condition_cov",h4("Coverage"),
       #           min=0.8,max=1,step=0.01,value = 0.9),
        br(),
        actionButton("run_button","Run Analysis",icon=icon("play")),
        
        h4("ReadMe"),
        h5("Identity: Two sequences are identified over 90% the same would be assumed as one cluster"),
        h5("Coverage: Two seqences length are covered over 90% each other"),
        img(src="coverage.png",height=140,width=300,align="left"),
        #img(src="coverage1.png",height=1250,width=620,align="left"),
        br(),
	downloadButton('downloadData','Download')
      ),
      
      mainPanel(
        tabsetPanel(
        #tabPanel('ReadMe'),
        # tabPanel('Rawdata_Check',
        #          DT::dataTableOutput('rawcheck'),
        #          plotOutput('distPlot')),
        tabPanel('Result',
                 DT::dataTableOutput('contents'),
                 plotOutput('distPlot')),
        tabPanel('Plot',fluidRow(column(width=12,plotOutput("plot"))),
                 fluidRow(column(width=12,plotOutput("plot1"))),
                 fluidRow(column(width=12,plotOutput("plot2")))
                                 
        )
        )
      )
      )
  ),
  tabPanel(
    "Test"
  )
  )
)

)

#source("/Users/tina/Desktop/CTMbio/Needs/MMseq2/MMseq2/script_final/cluster_function.R")
source("/home/shiny/script/modify_Abnum.R")
source("/home/shiny/script/runcluster.R")

server <- function(input, output, session){
  #set save dir
  savedir<-"/data/shiny/abcluster"
  rv <- reactiveValues(data = NULL)
  # data input
  myData <- eventReactive(input$file1,{
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    if (stringr::str_ends(inFile$datapath,'csv')) {
    data<-read.csv(inFile$datapath, header = TRUE)
    #data$newID<-paste0("ID_",c(1:dim(data)[1]))	
    rv$data<-data
    }else if (stringr::str_ends(inFile$datapath,"(xlsx|xls)")) {
    data<-read.xlsx(inFile$datapath)
    #data$newID<-paste0("ID_",c(1:dim(data)[1]))
    rv$data<-data
    }else if (stringr::str_ends(inFile$datapath,"txt")){
    data<-read.table(inFile$datapath, header=TRUE, sep="\t")
    #data$newID<-paste0("ID_",c(1:dim(data)[1]))
    rv$data<-data
}
  })
  
  newdir <- reactiveValues(dirpalce=NULL,fasta=data.frame(),dupset=data.frame(),
                           aa=data.frame(),outdata=data.frame(),all_out=data.frame(),all_out2=data.frame())
  
  # button invoke
  mmseqid<-eventReactive(input$run_button,{isolate(input$condition_id)})
  #mmseqcov<-eventReactive(input$run_button,{isolate(input$condition_cov)})
    
  #fr1<-"EIVLTQSPGTLSLSPGERATLSC"
  #fr2<-"IGWFRQAPGKEREGVAW"
  #fr3<-"YYADSVKGRFTISSDNAKNTVSLQMNSLKPEDTAKYYC"

  #Heavy Chain Framework seq (A0A0C4DH42 (HV366_HUMAN))
  #hfr1<-"EVQLVESGGGLIQPGGSLRLSCAASGFTVS"
  #hfr2<-"WVRQAPGKGLEWVS"
  #hfr3<-"RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
  #Light Chain Framework seq (>pdb|3FZU|L Chain L, immunoglobulin IgG1 Fab, light chain)
  #lfr1<-"DIQMTQSPSSLSASVGDRVTITC"
  #lfr2<-"WYQQKPGKAPKLLIY"
  #lfr3<-"GVPSRFSGSGSGTEFTLTVSSLQPEDFATYYC"  

  #save file
  
   observeEvent(input$run_button,{
    showModal(modalDialog(h2("Running..."),footer=NULL))
    #create dir
    inFile<-input$file1
    getname<-gsub("@.*","",inFile$name)
    getdirname<-gsub("\\..*","",gsub(".*@","",inFile$name))
    timeline<-gsub(" ","_",system("time date", TRUE))
    newdir$dirpalce<-paste0(savedir,"/",getname,"/",getdirname,"_",timeline)
    dir.create(newdir$dirpalce)
    # 
    # setwd(newdir)
    write.table(myData(),file=paste0(newdir$dirpalce,"/saveout.txt"),row.names = F,quote = F, sep="\t")
    #mmseqs2_analysis(myData(),mmseqid(),mmseqcov(),0,1,newdir)
    #mmseqs2_analysis(myData(),1,newdir$dirpalce)
    system(paste0("cd ",newdir$dirpalce),TRUE)
    
    if(input$checkGroup %in% c("vh","vhvlcdr")){
     if(input$checkGroup=="vh"){ 
      #1. =====HC only
     runcluster(myData(),newdir$dirpalce,mmseqid(),1,mode="vh")	
     outdata<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_Output.xlsx"),sheet=1)
     outdata$nH_CDR1<-nchar(outdata$H_CDR1)
     outdata$nH_CDR2<-nchar(outdata$H_CDR2)
     outdata$nH_CDR3<-nchar(outdata$H_CDR3)
     newdir$all_out<-unique(outdata)
     } 
    if (input$checkGroup=="vhvlcdr"){
    
    #2. =====HC_LC
    runcluster(myData(),newdir$dirpalce,mmseqid(),1,mode="vhvlcdr") 
    vhvl<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_HCLC_Result.xlsx"))
    vh<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_HC_Result.xlsx"))
    vl<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_LC_Result.xlsx"))
  
    vh<-vh[,c("ID","H_cluster")]
    vl<-vl[,c("ID","L_cluster")]
    
    all_out<-merge(vhvl,vh,by=c("ID"),all.x=T)
    all_out<-merge(all_out,vl,by=c("ID"),all.x=T)
    all_out<-unique(all_out)
    all_out$nH_CDR1<-nchar(all_out$H_CDR1)
    all_out$nH_CDR2<-nchar(all_out$H_CDR2)
    all_out$nH_CDR3<-nchar(all_out$H_CDR3)
    all_out$nL_CDR1<-nchar(all_out$L_CDR1)
    all_out$nL_CDR2<-nchar(all_out$L_CDR2)
    all_out$nL_CDR3<-nchar(all_out$L_CDR3)
    newdir$all_out<-all_out 
   # write.xlsx(newdir$all_out,paste0(newdir$dirpalce,"/Cluster_Output.xlsx"))
    Outlist<-list()
    Outlist[[1]]<-all_out
    if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_HC.txt"))){
    duphc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_HC.txt"),sep="\t",header=T)
    Outlist[[2]]<-duphc    
    }else{
    duphc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[2]]<-duphc
    }  
   if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_HCLC.txt"))){
    duphclc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_HCLC.txt"),sep="\t",header=T)
    Outlist[[3]]<-duphclc
    }else{
    duphclc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[3]]<-duphclc
    }
   if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_LC.txt"))){
    duplc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_LC.txt"),sep="\t",header=T)
    Outlist[[4]]<-duplc
    }else{
    duplc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[4]]<-duplc
    }
    if(length(is.na(all_out$H_cluster))!=0){
      hlackid<-all_out$ID[which(is.na(all_out$H_cluster))]
      for(hid in hlackid){
        all_out$H_cluster[which(all_out$ID==hid)]<-duphc$H_cluster[which(duphc$ID==hid)]
        }
    }
    if(length(is.na(all_out$L_cluster))!=0){
      llackid<-all_out$ID[which(is.na(all_out$L_cluster))]
      for(lid in llackid){
        all_out$L_cluster[which(all_out$ID==lid)]<-duplc$L_cluster[which(duplc$ID==lid)]
      }
    }   
    
   Outlist[[1]]<-all_out
   names(Outlist)<-c("Cluster_All_Result","Duplicated_set_HC","Duplicated_set_HCLC","Duplicated_set_LC")
   newdir$all_out<-Outlist[[1]]
   newdir$all_out2<-Outlist
   write.xlsx(newdir$all_out2,paste0(newdir$dirpalce,"/Cluster_Output.xlsx"))

}
      
    }else if (input$checkGroup=="vhvl") {
      
    AnnotateCDR(myData(),"/home/shiny/script",newdir$dirpalce)
    datainput<-read.table(paste0(newdir$dirpalce,"/data_output.txt"),header=T,sep="\t")
    runcluster(datainput,newdir$dirpalce,mmseqid(),1,mode="vhvlcdr")
    vhvl<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_HCLC_Result.xlsx"))
    vh<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_HC_Result.xlsx"))
    vl<-read.xlsx(paste0(newdir$dirpalce,"/Cluster_LC_Result.xlsx"))

    vh<-vh[,c("ID","H_cluster")]
    vl<-vl[,c("ID","L_cluster")]

    all_out<-merge(vhvl,vh,by=c("ID"),all.x=T)
    all_out<-merge(all_out,vl,by=c("ID"),all.x=T)
    all_out<-unique(all_out)
    newdir$all_out<-all_out
    all_out$nH_CDR1<-nchar(all_out$H_CDR1)
    all_out$nH_CDR2<-nchar(all_out$H_CDR2)
    all_out$nH_CDR3<-nchar(all_out$H_CDR3)
    all_out$nL_CDR1<-nchar(all_out$L_CDR1)
    all_out$nL_CDR2<-nchar(all_out$L_CDR2)
    all_out$nL_CDR3<-nchar(all_out$L_CDR3)
   # write.xlsx(newdir$all_out,paste0(newdir$dirpalce,"/Cluster_Output.xlsx"))
    Outlist<-list()
    Outlist[[1]]<-all_out
    if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_HC.txt"))){
    duphc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_HC.txt"),sep="\t",header=T)
    Outlist[[2]]<-duphc
    }else{
    duphc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[2]]<-duphc
    }
   if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_HCLC.txt"))){
    duphclc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_HCLC.txt"),sep="\t",header=T)
    Outlist[[3]]<-duphclc
    }else{
    duphclc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[3]]<-duphclc
    }
   if (file.exists(paste0(newdir$dirpalce,"/duplicated_seqs_LC.txt"))){
    duplc<-read.table(paste0(newdir$dirpalce,"/duplicated_seqs_LC.txt"),sep="\t",header=T)
    Outlist[[4]]<-duplc
    }else{
    duplc<-data.frame(Noduplicated=c("Noduplicated"))
    Outlist[[4]]<-duplc
    }
    
    if(length(is.na(all_out$H_cluster))!=0){
    hlackid<-all_out$ID[which(is.na(all_out$H_cluster))]
    for(hid in hlackid){
      all_out$H_cluster[which(all_out$ID==hid)]<-duphc$H_cluster[which(duphc$ID==hid)]
    }
  }
  if(length(is.na(all_out$L_cluster))!=0){
    llackid<-all_out$ID[which(is.na(all_out$L_cluster))]
    for(lid in llackid){
      all_out$L_cluster[which(all_out$ID==lid)]<-duplc$L_cluster[which(duplc$ID==lid)]
    }
  }

   Outlist[[1]]<-all_out
   names(Outlist)<-c("Cluster_All_Result","Duplicated_set_HC","Duplicated_set_HCLC","Duplicated_set_LC")
   newdir$all_out<-Outlist[[1]]
   newdir$all_out2<-Outlist
   write.xlsx(newdir$all_out2,paste0(newdir$dirpalce,"/Cluster_Output.xlsx"))  
      
    }
    
    removeModal()  

    })
      # ######plot
  plotInput<-reactive({
      
      ddata<-newdir$all_out
     # ddata$ID<-paste0("A",ddata$ID)
     if(length(which(is.na(ddata$H_cluster)))!=0){ 
    ddata<-ddata[-which(is.na(ddata$H_cluster)),]
}else{
    ddata<-ddata
	}	
      d1 <- data.frame(from="origin", to=unique(ddata$H_cluster))
      d2 <- data.frame(from=ddata$H_cluster,to=ddata$ID)
      edges <- rbind(d1, d2)
      name <- unique(c(as.character(edges$from), as.character(edges$to)))
      vertices <- data.frame(
        name=name,
        group=c( rep(NA,(length(unique(ddata$H_cluster))+1)) , ddata$H_cluster)
      )
      # Create a graph object
      mygraph <- graph_from_data_frame(edges, vertices=vertices)

      p<-ggraph(mygraph, layout = 'dendrogram') +
        geom_edge_diagonal() +
        geom_node_text(aes( label=name, filter=leaf, color=group) , angle=90 , hjust=1, nudge_y=-0.05,parse=TRUE) +
        ylim(-.6, NA) +ggtitle("H_cluster")+theme_void()+coord_fixed(8)+
        theme(legend.position="none",plot.title = element_text(color="Black", size=12, face="bold.italic"))
      
  })
  
  plotInput1<-reactive({
    if(input$checkGroup!="vh"){
      ddata<-newdir$all_out

       if(length(which(is.na(ddata$L_cluster)))!=0){
    	ddata<-ddata[-which(is.na(ddata$L_cluster)),]
	}else{
   	ddata<-ddata
     	}
     # ddata$ID<-paste0("A",ddata$ID)
      d1 <- data.frame(from="origin", to=unique(ddata$L_cluster))
      d2 <- data.frame(from=ddata$L_cluster,to=ddata$ID)
      edges <- rbind(d1, d2)
      name <- unique(c(as.character(edges$from), as.character(edges$to)))
      vertices <- data.frame(
        name=name,
        group=c( rep(NA,(length(unique(ddata$L_cluster))+1)) , ddata$L_cluster)
      )
      # Create a graph object
      mygraph <- graph_from_data_frame(edges, vertices=vertices)
      
      p<-ggraph(mygraph, layout = 'dendrogram') +
        geom_edge_diagonal() +
        geom_node_text(aes( label=name, filter=leaf, color=group) , angle=90 , hjust=1, nudge_y=-0.05,parse=TRUE) +
        ylim(-.6, NA) +ggtitle("L_cluster")+theme_void()+coord_fixed(8)+
        theme(legend.position="none",plot.title = element_text(color="Black", size=12, face="bold.italic")
    
)}
  })
  
  plotInput2<-reactive({
    if(input$checkGroup!="vh"){
      ddata<-newdir$all_out
     # ddata$ID<-paste0("A",ddata$ID)
      d1 <- data.frame(from="origin", to=unique(ddata$HC_LC_cluster))
      d2 <- data.frame(from=ddata$HC_LC_cluster,to=ddata$ID)
      edges <- rbind(d1, d2)
      name <- unique(c(as.character(edges$from), as.character(edges$to)))
      vertices <- data.frame(
        name=name,
        group=c( rep(NA,(length(unique(ddata$HC_LC_cluster))+1)) , ddata$HC_LC_cluster)
      )
      # Create a graph object
      mygraph <- graph_from_data_frame(edges, vertices=vertices)
      
      p<-ggraph(mygraph, layout = 'dendrogram') +
        geom_edge_diagonal() +
        geom_node_text(aes( label=name, filter=leaf, color=group) , angle=90 , hjust=1, nudge_y=-0.05,parse=TRUE) +
        ylim(-.6, NA) + ggtitle("H_L_cluster") + theme_void()+coord_fixed(8)+
        theme(legend.position="none",plot.title = element_text(color="Black", size=12, face="bold.italic"))}
  })
  
  
      #ggsave(file=paste0(newdir$dirpalce,"/_C",length(unique(ddata$cluster)),".pdf"),limitsize = FALSE,width=pretty(round(length(unique(ddata$cluster))))[2]/2.5,height =6)
  observeEvent(input$run_button,{
  if(input$checkGroup!="vh"){
  
  ggsave(paste0(newdir$dirpalce,"/HC_cluster_plot.pdf"),plotInput(),dpi = 600,height=8,width=24)
  ggsave(paste0(newdir$dirpalce,"/HCLC_cluster_plot.pdf"),plotInput2(),dpi = 600,height=8,width=24)
  ggsave(paste0(newdir$dirpalce,"/LC_cluster_plot.pdf"),plotInput1(),dpi = 600,height=8,width=24)
  }else{
 ggsave(paste0(newdir$dirpalce,"/HC_cluster_plot.pdf"),plotInput(),dpi = 600,height=8,width=24)
  }

})

  output$plot<-renderPlot(bg="white",{
    print(plotInput())})
  output$plot1<-renderPlot({
    print(plotInput1())
  })
  output$plot2<-renderPlot({
    print(plotInput2())
  })
  output$contents <- DT::renderDataTable({
    DT::datatable(newdir$all_out[c(1:20),])
  })
  
  output$downloadData<-downloadHandler(
  filename= function(){
  paste0(basename(as.character(newdir$dirpalce)),".zip")
},
  content<-function(file){
  tmpdir<-tempdir()
  setwd(tempdir())
 
 if (input$checkGroup=="vh"){
   file.rename(paste0(newdir$dirpalce,"/Cluster_Output.xlsx"),
   paste0(newdir$dirpalce,"/Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"))

   file.copy(paste0(newdir$dirpalce,"/Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"),file)
   file.copy(paste0(newdir$dirpalce,"/HC_cluster_plot.pdf"),file)
   fs<-c(paste0("Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"),"HC_cluster_plot.pdf")
 }else{
  file.rename(paste0(newdir$dirpalce,"/Cluster_Output.xlsx"),
   paste0(newdir$dirpalce,"/Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"))
 
  file.copy(paste0(newdir$dirpalce,"/Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"),file)
  file.copy(paste0(newdir$dirpalce,"/HC_cluster_plot.pdf"),file)
  file.copy(paste0(newdir$dirpalce,"/HCLC_cluster_plot.pdf"),file)
  file.copy(paste0(newdir$dirpalce,"/LC_cluster_plot.pdf"),file)  
  fs<-c(paste0("Cluster_Output_id_",mmseqid(),"_cov_","1","_.xlsx"),"HC_cluster_plot.pdf","HCLC_cluster_plot.pdf","LC_cluster_plot.pdf")
} 
 zip::zipr(zipfile=file,files=paste0(newdir$dirpalce,"/",fs))

},

 contentType="application/zip"

)   
}

shinyApp(ui,server)
