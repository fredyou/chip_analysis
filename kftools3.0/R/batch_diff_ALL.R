batch_diff_ALL <-
function(datafile=datafile,list=list,combination=combination,chipType="ge",picName="picName",pairedFC,setScriptPath="path"){
# Load the needed R programs from your path
#if system is linux

print("Sourcing files start...");

if(!missing(setScriptPath)){
    sourcefiles=list.files(path=setScriptPath,pattern="*.r",ignore.case=T,full.names=T)
    if(length(sourcefiles)!=0 && sourcefiles!=""){
      sapply(sourcefiles,source,.GlobalEnv)
    }else{stop("No r files in the set Path!!")}
# }else{
#   if(grepl("linux",Sys.info()['sysname'],ignore.case=T)){
#   	source('/home/corona/gaoy/programs/work/stat_Flag.r')
#   	source('/home/corona/gaoy/programs/work/Filter_Flag.r')
#   	source('/home/corona/gaoy/programs/work/BatchFoldchange.r')
#   	source('/home/corona/gaoy/programs/work/BatchFoldchange_HTA.r')
#   	source('/home/corona/gaoy/programs/work/diff_up_down_gene_count_hta.r')
#   	source('/home/corona/gaoy/programs/work/PCA_plot.r')
#    	source('/home/corona/gaoy/programs/work/hclust_plot.r')
#   	
#   }else if(grepl("windows",Sys.info()['sysname'],ignore.case=T)){
#   	source('D:/Expression_Scripts/work/stat_Flag.r')
# 	source('D:/Expression_Scripts/work/Filter_Flag.r')
# 	source('D:/Expression_Scripts/work/BatchFoldchange.r')
# 	source('D:/Expression_Scripts/work/BatchFoldchange_HTA.r')
# 	source('D:/Expression_Scripts/work/diff_up_down_gene_count_hta.r')
# 	source('D:/Expression_Scripts/work/PCA_plot.r')
# 	source('D:/Expression_Scripts/work/hclust_plot.r')
# }else{
#     writeLines("Something is wrong with your system type.\nIf your system is neither linux nor windows, please comment this sentence")
#     stop();
#   }
}
print("Sourcing files End...")

#Global var

suffixes="_NS";
colname2="name";#name

colname1 ="";#ProbeName
suffixes2="";#flags

#assign colname1, suffixes2
if(grepl("ge",chipType,ignore.case=T)|| grepl("lnc.*|linc.*",chipType,ignore.case=T)){

	colname1 ="ProbeName";#ProbeName
	suffixes2="_flags";#flags

}else if(grepl(paste("affy","|","hta",sep=""),chipType,ignore.case=T)){
    
	colname1 ="Probe_Set_ID";#ProbeName
	suffixes2="_call";#flags

}else if(grepl("miRNA",chipType,ignore.case=T)){

	colname1 ="SystematicName";#ProbeName
	suffixes2="_flags";#flags

}else{stop("You ve assigned a wrong chipType!!ChipType could be ge,affy,hta,miRNA")}

#picName is the pic name of PCA, hcluster
#only PCA and hcluster function need this var
if(!missing(picName)){
    picName="Book1.txt_AP_selection"    
}else{if(picName==""){stop("If u assign a picName, it mustn`t be blank!!")}}

if(missing(list)){
  list=read.delim("list.tsv",colClasses=rep("character",3))
}else if(class(list)=="character"){
  list=read.delim(list)
}else if(class(list)=="data.frame"){
  list=list
}else{stop("list should be filename or dataframe!!")}

if(missing(combination)){
  combination=read.delim("Combination.txt",colClasses=rep("character",3))
}else if(class(combination)=="character"){
  combination=read.delim(combination)
}else if(class(combination)=="data.frame"){
  combination=combination
}else{stop("combination should be filename or dataframe!!")}

if(missing(datafile) && !grepl("lnc.*|linc.*",chipType,ignore.case=T) ){
  datafile=read.delim("Book1.txt")
}

if(missing(pairedFC)){pairedFC=FALSE}
if(length(unique(list$Group))<dim(list)[1]){
  print("From list I can see there are multi samples in a group")
  print("Let`s do paired foldchange calculation!")
  pairedFC=TRUE
}else{pairedFC=FALSE}

if(!grepl("hta",chipType,ignore.case=T)){
  #if not hta
  print(colname1);
  print(colname2);
  print(suffixes2);
  diffAllSub=function(datafile=datafile){
        if(pairedFC==FALSE){
          print("UnpairedFC!")
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          Paired = TRUE,
                          is.log = 1,
                          is.FC = 0
          )
          
          #  Filter by FC = 2
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '_fc2',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          Paired = TRUE,
                          is.log = 1,
                          is.FC = 1,
                          envelope = 1,
                          FC = 2
          )
          
          #  Filter by FC = 3
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '_fc3',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          Paired = TRUE,
                          is.log = 1,
                          is.FC = 1,
                          envelope = 1,
                          FC = 3
          )
        }else{#pairedFC
          print("use pairedFC")
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          is.log = 1,
                          is.FC = 0
          )
          
          #  Filter by FC = 2
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '_fc2',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          is.log = 1,
                          is.FC = 1,
                          envelope = 1,
                          FC = 2
          )
          
          #  Filter by FC = 3
          BatchFoldchange(DataFile = datafile,
                          ContentFile = list,
                          BatchInfoFile = combination,
                          OutfileSuffixes = '_fc3',
                          FilePath = './',
                          colname1 = colname1,
                          colname2 = colname2,
                          prefixes = '',
                          suffixes = '_NS',
                          suffixes2 = suffixes2,
                          is.log = 1,
                          is.FC = 1,
                          envelope = 1,
                          FC = 3
          )
    }	
    
  }
  
  if(grepl("lnc.*|linc.*",chipType,ignore.case=T)){
    if(exists("dataG")){ rm(dataG,envir=.GlobalEnv)}
    if(length(Sys.glob("Book*mRNA*.txt"))!=0 && length(Sys.glob("Book*lnc*.txt"))!=0){
      dir.create("mRNA",showWarnings=F)
      dir.create("lncRNA",showWarnings=F)
      #file.rename(from="Book_mRNA.txt",to="mRNA/Book_mRNA.txt")
      #file.rename(from="Book_lncRNA.txt",to="lncRNA/Book_lncRNA.txt")
      currentDir=getwd()
      setwd(paste(currentDir,"/mRNA",sep=""))
      diffAllSub("../Book_mRNA.txt")
      setwd(paste(currentDir,"/lncRNA",sep=""))
      if(exists("dataG")){ rm(dataG,envir=.GlobalEnv)}
      diffAllSub("../Book_lncRNA.txt")
      setwd(currentDir)
      
    }else{
      stop("If you r analysing lncRNA, you should have both Book_mRNA.txt and Book_lnc.txt right here!") 
    }
  }else{diffAllSub(datafile=datafile)}
  
 
  
    
}else{
  #else if chipType is hta
  if(pairedFC==FALSE){
  	BatchFoldchangeHTA(DataFile = datafile,
                ContentFile = list,
                BatchInfoFile = combination,
                OutfileSuffixes = '',
                FilePath = './',
                colname1 = colname1,
                colname2 = colname2,
                prefixes = '',
                suffixes = '_NS',
                suffixes2 = suffixes2,
                Paired = TRUE,
                is.log = 1,
                is.FC = 0
            )

	#  Filter by FC = 2
	BatchFoldchangeHTA(DataFile = datafile,
	                ContentFile = list,
	                BatchInfoFile = combination,
	                OutfileSuffixes = '_fc2',
	                FilePath = './',
	                colname1 = colname1,
	                colname2 = colname2,
	                prefixes = '',
	                suffixes = '_NS',
	                suffixes2 = suffixes2,
	                Paired = TRUE,
	                is.log = 1,
	                is.FC = 1,
	                envelope = 1,
	                FC = 2
	               )

	#  Filter by FC = 3
	BatchFoldchangeHTA(DataFile = datafile,
	                ContentFile = list,
	                BatchInfoFile = combination,
	                OutfileSuffixes = '_fc3',
	                FilePath = './',
	                colname1 = colname1,
	                colname2 = colname2,
	                prefixes = '',
	                suffixes = '_NS',
	                suffixes2 = suffixes2,
	                Paired = TRUE,
	                is.log = 1,
	                is.FC = 1,
	                envelope = 1,
	                FC = 3
	               )

  }else{
  	BatchFoldchangeHTA(DataFile = datafile,
                ContentFile = list,
                BatchInfoFile = combination,
                OutfileSuffixes = '',
                FilePath = './',
                colname1 = colname1,
                colname2 = colname2,
                prefixes = '',
                suffixes = '_NS',
                suffixes2 = suffixes2,
                is.log = 1,
                is.FC = 0
               )
                           
	#  Filter by FC = 2
	BatchFoldchangeHTA(DataFile = datafile,
	                ContentFile = list,
	                BatchInfoFile = combination,
	                OutfileSuffixes = '_fc2',
	                FilePath = './',
	                colname1 = colname1,
	                colname2 = colname2,
	                prefixes = '',
	                suffixes = '_NS',
	                suffixes2 = suffixes2,
	                is.log = 1,
	                is.FC = 1,
	                envelope = 1,
	                FC = 2
	               )

	#  Filter by FC = 3
	BatchFoldchangeHTA(DataFile = datafile,
	                ContentFile = list,
	                BatchInfoFile = combination,
	                OutfileSuffixes = '_fc3',
	                FilePath = './',
	                colname1 = colname1,
	                colname2 = colname2,
	                prefixes = '',
	                suffixes = '_NS',
	                suffixes2 = suffixes2,
	                is.log = 1,
	                is.FC = 1,
	                envelope = 1,
	                FC = 3
	               )

  }#end of else pairedFC
  dirs=dir(pattern="VS")[file.info(dir(pattern="VS"))$isdir]
  for(i in 1:length(dirs)){
    files=list.files(path=dirs[i],pattern="fc")
    diff_up_down_gene_count(datafile=files,filepathin=paste(dirs[i],"/",sep=""),filepathout="./")
  }
  scatterplot_NoFlag(combination = combination, list = list, datafile = datafile)
  PCA_plot(datafile=datafile,contentfile = list,suffixes2='_NS',filepathin="./",filepathout="./")
  hclust_plot(datafile=datafile,contentfile = list,suffixes2='_NS',filepathin="./",filepathout="./")
}#end of else hta

}
