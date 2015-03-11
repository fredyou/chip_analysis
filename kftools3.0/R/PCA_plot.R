PCA_plot <-
function(datafile="Book1.txt",contentfile = 'list.tsv',suffixes2='_NS',filepathin="./",filepathout="./",picName="g1_vs_g2"){

group_tb=""
if(class(contentfile)=="character"){
group_tb <-  read.delim(paste(filepathin,contentfile,sep=''),sep="\t",quote="",header=T,fill=T,stringsAsFactors=F); 
}else{group_tb=contentfile;contentfile=paste(contentfile[1,3],"_VS_",contentfile[2,3]);}
tb1=""
if(class(datafile)=="character"){
tb1 <- read.delim(paste(filepathin,datafile,sep=''),sep="\t",quote="",header=T,fill=T,stringsAsFactors=F);
}else{
tb1=datafile;datafile=contentfile;
}
# rownames(tb1) <- tb1[,"ProbeName"]
tb2 <- tb1[,grep(suffixes2,colnames(tb1))]
colnames(tb2) <- gsub(suffixes2,"",colnames(tb2))
## PCA
fit.pca <- princomp(tb2)
var_prop <- fit.pca$sdev^2/sum(fit.pca$sdev^2) 
sam_names    <- rownames(fit.pca$loadings)
sam_colors  <- as.factor(group_tb[match(rownames(fit.pca$loadings), group_tb[,"name"]),"Group"])
pch_pca <- as.character(sam_colors)
for(i in unique(pch_pca)){
pch_pca[!is.na(match(pch_pca,i))] <- 1:length(which(!is.na(match(pch_pca,i))))
} 
pch_pca <- as.integer(pch_pca) 

if(!missing(picName)){
    datafile=picName
}else{
  if(class(datafile)!="character"){datafile="Book1.txt"}
}
## sam_colors <- 1:length(sam_names) 
pdf(file=paste(filepathout,datafile,"_PCA_plot_2D.pdf",sep=""))
par(xpd=NA,oma=c(0,0,0,5)) 
plot(fit.pca$loadings[,1:2]
              ,xlab=paste("Comp 1", round(var_prop[1],3),sep=" : ")
              ,ylab=paste("Comp 2", round(var_prop[2],3),sep=" : ")
              ,col=sam_colors
              ,cex=1
              #,box.col="green",  fill="red"
              ,pch=pch_pca
              ,main=gsub(".txt","",datafile))              
## text(fit.pca$loadings[,1:2], rownames(fit.pca$loadings),col=as.numeric(as.factor(group_tb[match(rownames(fit.pca$loadings), group_tb[,"name"]),"Group"])),cex=1.0)  
## legend(max(fit.pca$loadings[,1]),max(fit.pca$loadings[,2]),sam_names,sam_colors) 
## sam_names[1] <- paste(rep("a",15),collapse="")
legend(par("usr")[2] + (max(fit.pca$loadings[,1])-min(fit.pca$loadings[,1]))/50   ,par("usr")[4],legend=sam_names
        ,col=sam_colors,pch=pch_pca,xjust=0, yjust=1.0,ncol=1,cex=min(1,8/max(nchar(sam_names))))                               


dev.off()

#### png output
png(file=paste(filepathout,datafile,"_PCA_plot_2D.png",sep=""),height=14,width=14,units = "cm",res=600)
par(xpd=NA,oma=c(0,0,0,5)) 
plot(fit.pca$loadings[,1:2]
              ,xlab=paste("Comp 1", round(var_prop[1],3),sep=" : ")
              ,ylab=paste("Comp 2", round(var_prop[2],3),sep=" : ")
              ,col=sam_colors
              ## ,col=1:length(sam_names)
              ,cex=1
              #,box.col="green",  fill="red"
              ,pch=pch_pca
              ,main=gsub(".txt","",datafile))              
## text(fit.pca$loadings[,1:2], rownames(fit.pca$loadings),col=as.numeric(as.factor(group_tb[match(rownames(fit.pca$loadings), group_tb[,"name"]),"Group"])),cex=1.0)  
## legend(max(fit.pca$loadings[,1]),max(fit.pca$loadings[,2]),sam_names,sam_colors) 
## sam_names[1] <- paste(rep("a",15),collapse="")
legend(par("usr")[2] + (max(fit.pca$loadings[,1])-min(fit.pca$loadings[,1]))/50   ,par("usr")[4],legend=sam_names
        ,col=sam_colors,pch=pch_pca,xjust=0, yjust=1.0,ncol=1,cex=min(1,8/max(nchar(sam_names))))                               

dev.off()

   }
