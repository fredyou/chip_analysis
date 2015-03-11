Batch_ttestHTA <-
function(datafile='Book1.txt',
                      contentfile = 'list.tsv',
                      Batch_infofile = 'Combination.txt',
                      paired_ttest_fname= 'paired.txt',
                      filepathin= '',
                      filepathout = '',
                      colname1 = 'ProbeName',
                      colname2 = 'filename',
                      prefixes = '',
                      suffixes = '_NS',
                      
                      T.type = 2,
                      is.log = 1,
                      cutoff = 1 ,
                      is.adjustP = 0)
{
  #-------------Read data--------------#
  #read expression file

  data="";
  if(class(datafile)=="character"){
            data<-read.delim(paste(filepathin,datafile,sep=''),sep="\t",quote="\"",header=T,fill=T);
            }else{
                      data=datafile
                                rm(datafile)    
  }  
    
    
  probeset<-as.vector(data[[colname1]])  
  
  #read Batch information file
  #Content of Combination.txt=Conditions
  Conditions="";
    if(class(Batch_infofile)=="data.frame"){Conditions=Batch_infofile}else{
                 Conditions <- read.delim(file.path(filepathin, Batch_infofile), sep = "\t", header = T, fill = T,check.names=FALSE)
                                     
                       } 
      
      
  Batch_info<-subset(Conditions,Selected=='Y')
  #print 
  print(paste("The number of iterations is ",as.character(nrow(Batch_info))))
  
  #read sample file
  sampleinfo=""
    if(class(contentfile)=="data.frame"){ sampleinfo=contentfile }else{
               sampleinfo <- read.delim(file.path(filepathin, contentfile), sep = "\t", header = T, fill = T,check.names=FALSE)
                        
    }

  #content of list.tsv
  
  data_annot<-data[,-grep(suffixes,colnames(data))]
  if (ncol(as.matrix(data_annot))==1){
    is.annot='N' 
  } else {
    is.annot='Y';
    data_annot<-as.matrix(data_annot[,-1]);
  }#
  
  for (i in 1:nrow(Batch_info)){
    
    #read label file
    #indicating elements or rows to keep: missing values are taken as false
    #Below is to subset sampleinfo by Group
    label1_raw <-subset(sampleinfo,Group==as.character(Batch_info[i,1]),select=colname2)
    label1_raw<-as.matrix(label1_raw)
    label1<-paste(prefixes,label1_raw,suffixes,sep="")
    
    label2_raw<-subset(sampleinfo,Group==as.character(Batch_info[i,2]),select=colname2)
    label2_raw<-as.matrix(label2_raw)
    label2<-paste(prefixes,label2_raw,suffixes,sep="")
    
    #read Flag
#     label1_Flag<-paste(prefixes,label1_raw,suffixes2,sep="");
#     label2_Flag<-paste(prefixes,label2_raw,suffixes2,sep="");
#     
    #print group
    print(paste("the NO is",as.character(i)))
    print(as.character(Batch_info[i,1]))
    print(label1_raw)
    print(as.character(Batch_info[i,2]))
    print(label2_raw)
    
    #------------t-test------------#
    data1<-subset(data,select=c(label1))
    data2<-subset(data,select=c(label2))
    ## if paired ttest
    ##
    result<-T_testHTA(data1,data2,T_type=T.type,islog=is.log,isadjustP=is.adjustP)		
    if (is.annot=='N')
      result<-cbind(probeset,result) else
        result<-cbind(probeset,result,data_annot);
    
    result <- subset(result,pvalues<=cutoff)
    colnames(result)[1]=colname1
    
    #-------------output result-----------------#
    #create output result directory with annotation
    ## outputdir<-paste(filepathout,Batch_info[i,1],"_VS_",Batch_info[i,2],sep='')
    ##dir.create(outputdir)
    outputdir <- filepathout
    assign(paste(Batch_info[i,1],"_VS_",Batch_info[i,2],sep=''),result,envir=.GlobalEnv)
    write.table(result,file=paste(outputdir,"/",Batch_info[i,1],"_VS_",Batch_info[i,2],".txt",sep=''),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  }
  
}
