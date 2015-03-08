prettyPCAPlot <-
function(datafile="Book1.txt",list="list.tsv",suffix="_NS",s3D="FALSE",type="percent",picName="Book1_selectAP"){
  library(gridExtra)
  library(gtable)
  library(grid)
  library(plyr)
  library(ggplot2)
  library(grDevices)
  pcaType=""
  if(!missing("type")){
    if(type=="percent"){ pcaType="percent"}else{pcaType="score"}
  }else{
    pcaType="percent"   
  }
  if(missing(picName)){
    if(class(datafile)=="character"){picName=datafile}
    else{picName="Book1.txt_AP_selection"}
  }else{
    picName=picName
  }
     
  #Processing Input Data
  if(class(list)=="character"){
    list=read.delim(list)
  }else{list=list}
  if(class(datafile)=="character"){
    data=read.delim(datafile,na.strings="",quote="",sep="\t",stringsAsFactors=F)

  }else{data=datafile}
  rownames(list)=paste(list$name,"_NS",sep="")
  numGroups=length(unique(list$Group))
  
  numProbe=grep("probe|systematic",colnames(data),ignore.case=T)[1]
  rownames(data)=data[,numProbe]
  #data=data[,grep(suffix,colnames(data))]
  data=data[,paste(list$name,suffix,sep="")]
  
  
  
  #Color hue and shape generate
  #colorhue with 128 distinctive colors
  #colorhue with 12 distinctive colors
  #In theory , 300 samples need 12 colors with 25 pch (12*25=300) to discriminate.
  #colorhue=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6as3D9a","#ffff99","#b15928")
  
      colorhue="";
      if(numGroups<=12){
        colorhue=c("#a6cee3",  "#1f78b4"	,"#b2df8a"	,"#33a02c"	,"#fb9a99"	,"#e31a1c"	,"#fdbf6f"	,"#ff7f00"	,"#cab2d6"	,"#6a3d9a"	,"#ffed6f","#b15928")
        #colorhue=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6as3D9a","#ffff99","#b15928")
      }else{
        #122 in total
        colorhue=c("#D063AD",  "#4DF036",  "#E4952E",	"#66DDDB",	"#232E34",	"#34711F",	"#4D85D9",	"#E8C1B5",	"#95241C",	"#D8E638",	"#DF49EF",	"#59E193",	"#EC2171",	"#6B4E1A",	"#521C4E",	"#75A279",	"#8153BF",	"#4D8BA4",	"#9E5F6A",	"#4B1819",	"#E4E39B",	"#E95A20",	"#CEB6E1",	"#42A82A",	"#A2A92F",	"#EE605A",	"#8C66A5",	"#E29F77",	"#D328AC",	"#E7C22E",	"#89EE73",	"#324210",	"#9D8028",	"#244C6B",	"#AF1E70",	"#3B6968",	"#64EEC5",	"#989F55",	"#D26FE9",	"#A1ABAB",	"#E07188", "#83273C",	"#574235","#3F6C42","#E093E4",	"#56A19A",	"#5DC5E3",	"#C9E268",	"#746D46",	"#666E1F",	"#693685",	"#5A515E",	"#64A7E3",	"#9565F0","#DAAF64",	"#E1322E",	"#A3A1EA",	"#65B75A",	"#302749","#C3869F",	"#7C8275",	"#C07D3B",	"#263C28",	"#7B4F76",	"#8B869E",	"#C76548",	"#9ABADB",	"#B084EB",	"#6780EE",	"#A77557",	"#B8B698",	"#637EB5",	"#9FE694",	"#4A9E61",	"#6C9B31",	"#4ABD97",	"#DDD2E1",	"#C25922",	"#E38275",	"#E1C65A",	"#DB558D",	"#AF2E5C",	"#DFDEC9",	"#9AF034",	"#3B182F",	"#BF3A49",	"#EC87B9",	"#A75379",	"#E1BF93",	"#E8A1A4",	"#EA4BAA",	"#3A2710",	"#796BBF",	"#EFB6D3",	"#864744",	"#2D2960",	"#B1D27D",	"#973F97",	"#A248C5",	"#A49C6B",	"#54D638",	"#E3A5DD",	"#652E49","#C0EAB5",	"#4C5B87",	"#AF807F",	"#DF6BCC",	"#44D763",	"#7B6862",	"#A872C4",	"#4C5DCA",	"#AC7BA9",	"#E73187",	"#A593C9",	"#E241CF",	"#AE309C",	"#BE9E9D",	"#62417D",	"#982F72",	"#73215C",	"#B564A4",	"#B72785")    
      }
      colors=sample(colorhue,numGroups)
      names(colors)=unique(as.character(list$Group))
      
      
      #initialize
      ncolused=numGroups
      #pchGroupNumber=""
      
      shapes=""
      countGroups=count(list,"Group")
      maxGroupCount=max(countGroups$freq)
     
      if(dim(list)[1]>23){
         #pchGroupNumber=floor(numGroups/ncolused)
        if(maxGroupCount>23){
          shapes= unlist(sapply(countGroups$freq,function(x){ if(x<=44){sample(c(0:20,23,25,33:47,60:64),x)}else{sample(c(0:20,23,25,33:126),x,replace=T)} }));
        }else{shapes=unlist(sapply(countGroups$freq,function(x){ sample(c(0:20,23,25,35,38),x) }))} 
        
      }else{
        shapes= sample(c(0:20,23,25),dim(list)[1]);
      }
     
     
      #For extention of this function, you can assign dataG to Globenv
      #Then in this function write: if(exists(dataG)){....}else{ for(f_name in data){}  }
  
  #Ploting begin
  if(pcaType=="score"){
    
      #for(f_name in datafile){
        #system( "dos2unix f_name" )
        #data=read.delim(f_name)
         
        rawdata=2^data
        
        m=t(as.matrix(rawdata))
        
        pca <- prcomp(m,scale=F)
        scores <- as.data.frame(pca$x[,1:4])#PC1-4 
        sam_scores=rownames(scores)
        Samples=as.factor(list[sam_scores,"name"])
         
        
        colorplatte=unname(colors[list[sam_scores,"Group"]])
        
        
        if(dim(list)[1]<=32){
         p=ggplot(data=scores,aes(x=PC1,y=PC2,colour=Samples,shape=Samples))+geom_point(size=3)+
           geom_vline(xintercept=0)+geom_hline(yintercept=0)+scale_shape_manual(values=shapes)+
           scale_colour_manual(values=colorplatte)+
           scale_x_continuous(breaks = round(seq(min(scores$PC1), max(scores$PC1), by = max(scores$PC1)/2),1))+
           theme_classic()
           #+coord_cartesian(xlim=c(min(scores$PC1)-5000,max(scores$PC1)+5000))
         ggsave(filename=paste(paste(unique(list$Group),collapse="_"),"_PCA.png",sep=""),plot=p,dpi=300)
         ggsave(filename=paste(paste(unique(list$Group),collapse="_"),"_PCA.pdf",sep=""),plot=p,dpi=300)
        }else{
          #If there are too many samples
          #Legends should be drawn seperatly
          #calculate how many lines or columns to draw a legend
          devheight=par("fin")[2] #inches
          devwidth=par("fin")[1] #inches
          numchars=ceiling(mean(nchar(as.character(list$name))))
          nlegCols=par("din")[1]/(numchars*par("csi")*0.35)      
          nlegCols=floor(nlegCols)      
          nlegRows=ceiling(dim(list)[1]/nlegCols)
          p=ggplot(data=scores,aes(x=PC1,y=PC2,colour=Samples,shape=Samples))+geom_point(size=3)+
            geom_vline(xintercept=0)+geom_hline(yintercept=0)+scale_shape_manual(values=shapes)+scale_colour_manual(values=colorplatte)+
            scale_x_continuous(breaks = round(seq(min(scores$PC1), max(scores$PC1), by = max(scores$PC1)/2),1))+theme(legend.position="bottom")+guides(colour=guide_legend(nrow = nlegRows,ncol=nlegCols),shape=guide_legend(nrow = nlegRows,ncol=nlegCols))+
            theme_classic()
          panel=gtable_filter(ggplot_gtable(ggplot_build(p)), "guide-box")
          leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p)), "guide-box") 
          leg1$heights=panel$heights
          print(paste("gtable legend height:",convertUnit(leg1$height,unitTo="inch"),sep=""))      
          #plot(leg1)
          p=ggplot(data=scores,aes(x=PC1,y=PC2,colour=Samples,shape=Samples))+geom_point(size=3)+
            geom_vline(xintercept=0)+geom_hline(yintercept=0)+scale_shape_manual(values=shapes)+
            scale_colour_manual(values=colorplatte)+
            scale_x_continuous(breaks = round(seq(min(scores$PC1), max(scores$PC1), by = max(scores$PC1)/2),1))+theme(legend.position="none")+guides(colour=guide_legend(nrow = nlegRows,ncol=nlegCols),shape=guide_legend(nrow = nlegRows,ncol=nlegCols))+
            theme_classic()
          #unit.c creates a list of unit object
          
    #       plotNew <- arrangeGrob(leg1, p, heights = unit.c(leg1$height, unit(1, "npc") - leg1$height))
    #       grid.newpage()
    #       grid.draw(plotNew)
    #        plotNew=arrangeGrob(p+theme(legend.position="none"),leg1, nrow=30)
    #       
          #Two pictures would have 2* par("din")[2]=185mm *2=370mm if been stacked, so png height should be no less than 370mm 
          png( paste(paste(unique(list$Group),collapse="_"),"_PCA.png",sep=""),width = 480, height = as.vector(convertUnit(panel$heights,"mm"))*2+100, units = "mm",res=300)
          grid.arrange(p,leg1,nrow=2)
          dev.off()
    #       grid.newpage()
    #       pushViewport(viewport(layout = grid.layout(1,2)))#x:2elements,y:2elements
    #       vplayout = function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}
          
          pdf(paste(paste(unique(list$Group),collapse="_"),"_PCA.pdf",sep=""),width=devwidth+3,height=devheight*2+3)
          grid.arrange(p,leg1,nrow=2)
          dev.off()
        }#end of else in pcaType score
      #}#end of for
  }#End of if pcaType is score
  else{
    require(gplots)
    library(fpc)
    library(rgl)
    #for(f_name in datafile){
      #system( "dos2unix f_name" )
      #data=read.delim(f_name)
       
      fit.pca <- princomp(data)
    
      var_prop <- fit.pca$sdev^2/sum(fit.pca$sdev^2)
      
      Samples=as.factor(list[rownames(fit.pca$loadings),"name"])
     #sam_fit=gsub("_NS","",rownames(fit.pca$loadings))
      cat("Samples: ",as.character(Samples))
      colorplatte=unname(colors[list[rownames(fit.pca$loadings),"Group"]])
      result=as.data.frame(fit.pca$loadings[,1:2],Samples=Samples)
      if(dim(list)[1]<=32){
        p=ggplot(data=result,aes(x=Comp.1,y=Comp.2,colour=Samples,shape=Samples))+
          geom_point(size=5)+scale_shape_manual(values=shapes)+
          scale_colour_manual(values=colorplatte)+
          theme_classic()+
          theme(panel.background = element_rect(colour = "black", size=2),plot.title=element_text(family="serif",face="bold",size=18))
        #+
        #  geom_text(aes(label=Samples,alpha=0.7),size=3,hjust=-0.1,vjust=-1,angle=45,show_guide=F)
        #+ggtitle(paste(paste(unique(list$Group),collapse="_"),"_PCA",sep=""))
        #+coord_cartesian(xlim=c(min(scores$PC1)-5000,max(scores$PC1)+5000))
        ggsave(filename=paste(paste(unique(list$Group),collapse="_"),"_PCA.png",sep=""),plot=p,dpi=300)
        ggsave(filename=paste(paste(unique(list$Group),collapse="_"),"_PCA.pdf",sep=""),plot=p,dpi=300)
        
      }else{
        devheight=par("fin")[2] #inches
        devwidth=par("fin")[1] #inches
        numchars=ceiling(mean(nchar(as.character(list$name))))
        nlegCols=par("din")[1]/(numchars*par("csi")*0.35)      
        nlegCols=floor(nlegCols)      
        nlegRows=ceiling(dim(list)[1]/nlegCols)
        p=ggplot(data=result,aes(x=Comp.1,y=Comp.2,colour=Samples,shape=Samples))+
          geom_point(size=5)+scale_shape_manual(values=shapes)+
          theme(legend.position="bottom")+
          scale_colour_manual(values=colorplatte)+
          guides(colour=guide_legend(nrow = nlegRows,ncol=nlegCols),shape=guide_legend(nrow = nlegRows,ncol=nlegCols))+
          theme_classic()
          #+coord_cartesian(xlim=c(min(scores$PC1)-5000,max(scores$PC1)+5000))
        panel=gtable_filter(ggplot_gtable(ggplot_build(p)), "guide-box")
        leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p)), "guide-box") 
        leg1$heights=panel$heights
        leg1$widths=panel$width
        print(paste("gtable legend height:",convertUnit(leg1$height,unitTo="inch"),sep=""))      
        #plot(leg1)
        p=ggplot(data=result,aes(x=Comp.1,y=Comp.2,colour=Samples,shape=Samples))+
          geom_point(size=5)+scale_shape_manual(values=shapes)+
          scale_colour_manual(values=colorplatte)+theme(legend.position="none",panel.background = element_rect(colour = "black", size=2),plot.title=element_text(family="serif",face="bold",size=18))+
          theme_classic()
        #+
        #  geom_text(aes(label=Samples,alpha=0.7),size=3,hjust=-0.1,vjust=-1,angle=45,show_guide=F)
        #+ggtitle(paste(paste(unique(list$Group),collapse="_"),"_PCA",sep=""))
        #unit.c creates a list of unit object
        
        #       plotNew <- arrangeGrob(leg1, p, heights = unit.c(leg1$height, unit(1, "npc") - leg1$height))
        #       grid.newpage()
        #       grid.draw(plotNew)
        #        plotNew=arrangeGrob(p+theme(legend.position="none"),leg1, nrow=30)
        #       
        #Two pictures would have 2* par("din")[2]=185mm *2=370mm if been stacked, so png height should be no less than 370mm 
        png( paste(paste(unique(list$Group),collapse="_"),"_PCA.png",sep=""),width =388, height = as.vector(convertUnit(panel$heights,"mm"))*2+100, units = "mm",res=300)
        grid.arrange(p,leg1,nrow=2)
        dev.off()
        #       grid.newpage()
        #       pushViewport(viewport(layout = grid.layout(1,2)))#x:2elements,y:2elements
        #       vplayout = function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}
        
        pdf(paste(paste(unique(list$Group),collapse="_"),"_PCA.pdf",sep=""),width=devwidth+3,height=devheight*2+3)
        grid.arrange(p,leg1,nrow=2)
        dev.off()
      }  
    #}#end of for loop in else  
  }#End of else pcaType is percent
}
