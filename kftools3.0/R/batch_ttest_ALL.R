batch_ttest_ALL <-
function(datafile=datafile,list=list,combination=combination,paired,chipType="chipType",picName="picName",setScriptPath="scriptPath"){
# Load the needed R programs from your path
#if system is linux
  
if(!missing(setScriptPath)){
    sourcefiles=list.files(path=setScriptPath,pattern="*.r",ignore.case=T,full.names=T)
    if(length(sourcefiles)!=0 && sourcefiles!=""){
      sapply(sourcefiles,source,.GlobalEnv)
    }else{stop("No r files in the set Path!!")}
}else{
  if(grepl("linux",Sys.info()['sysname'],ignore.case=T)){
  	source('/home/corona/gaoy/programs/work/T_test.r')
  	source('/home/corona/gaoy/programs/work/Batch_ttest.r')
  	source('/home/corona/gaoy/programs/work/stat_Flag.r')
  	source('/home/corona/gaoy/programs/work/Filter_Flag.r')
  	source('/home/corona/gaoy/programs/work/Batch_filter.r')
  	#source('/home/corona/gaoy/programs/work/Batch_volcanoplot.r')
  	source('/home/corona/gaoy/programs/work/prettyVolcanoPlot.R')
  	source('/home/corona/gaoy/programs/work/heatmap_plot.r')
  	source('/home/corona/gaoy/programs/work/prettyheatmap_plot.r')
  	source('/home/corona/gaoy/programs/work/PCA_plot.r')
  	#source('/home/corona/gaoy/programs/work/prettyPCAPlot.R')
    source('/home/corona/gaoy/programs/work/hclust_plot.r')
  	source('/home/corona/gaoy/programs/work/Total_AP_selection.r')
  	source('/home/corona/gaoy/programs/work/diff_up_down_gene_count.r')
  	#source('/home/corona/gaoy/programs/work/code.R')
  	source('/home/corona/gaoy/programs/work/cor_plot.r')
  	source('/home/corona/gaoy/programs/work/T_test_HTA.r')
  	source('/home/corona/gaoy/programs/work/Batch_ttest_HTA.r')
  	source('/home/corona/gaoy/programs/work/Batch_filter_HTA.r')
  }else if(grepl("windows",Sys.info()['sysname'],ignore.case=T)){
  	source('D://Expression_Scripts//work//T_test.r')
  	source('D://Expression_Scripts//work//Batch_ttest.r')
  	source('D://Expression_Scripts//work//stat_Flag.r')
  	source('D://Expression_Scripts//work//Filter_Flag.r')
  	source('D://Expression_Scripts//work//Batch_filter.r')
  	#source('D://Expression_Scripts//work//Batch_volcanoplot.r')
  	source('D://Expression_Scripts//work//prettyVolcanoPlot.R')
    source('D://Expression_Scripts//work//heatmap_plot.r')
  	source('/home/corona/gaoy/programs/work/prettyheatmap_plot.r')
  	source('D://Expression_Scripts//work//PCA_plot.r')
  	#source('D://Expression_Scripts//work//prettyPCAPlot.R')
    source('D://Expression_Scripts//work//hclust_plot.r')
  	source('D://Expression_Scripts//work//Total_AP_selection.r')
  	source('D://Expression_Scripts//work//diff_up_down_gene_count.r')
  	#source('D://Expression_Scripts//work//code.R')
  	source('D://Expression_Scripts//work//cor_plot.r')
  	source('D://Expression_Scripts//work//T_test_HTA.r')
  	source('D://Expression_Scripts//work//Batch_ttest_HTA.r')
  	source('D://Expression_Scripts//work//Batch_filter_HTA.r')
  }else{
    writeLines("Something is wrong with your system type.\nIf your system is neither linux nor windows, please comment this sentence")
    stop();
  }
}


#Global var
suffixes="_NS";
colname2="name";#name

colname1 ="";#ProbeName
suffixes2="";#flags

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

if(missing(datafile) && !grepl("lnc.*|linc.*",chipType,ignore.case=T)){
  datafile=read.delim("Book1.txt",na.strings="",quote="",sep="\t",stringsAsFactors=F)
}

#paired T or not ;1 is paired
ttype=1
if(!missing(paired)){ttype=2}

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
if(missing(picName)){
    picName="Book1.txt_AP_selection"    
}else{if(picName==""){stop("If u assign a picName, it mustn`t be blank!!")}}

if(!grepl("hta",chipType,ignore.case=T)){
  #if not hta
  ttestSub=function(datafile=datafile){
      Batch_ttest(datafile=datafile,
                  contentfile=list,
                  Batch_infofile=combination,
                  paired_ttest_fname= 'paired.txt',
                  filepathin="./",
                  filepathout="./",
                  colname1 = colname1,
                  colname2='name',
                  prefixes='',
                  suffixes=suffixes,
                  suffixes2=suffixes2,
                  T.type = ttype,
                  is.log = 1,
                  cutoff = 1,
                  is.adjustP = 0
      )
      
      # Filter by P-value < 0.05
      Batch_filter(
        Batch_infofile = combination,
        contentfile = list,
        outfile_suffixes = '_p005',
        filepath="./",
        colname2 = 'name',
        prefixes = '',
        suffixes1 = suffixes,
        suffixes2 = suffixes2,
        P_value = 0.05,
        is.FC = 0,
        FC = 2 ,
        envelope = 1,
        is.adjustP = 0,
        is.mean_sel=0,
        mean_sel=6
      )
      
      
      # Filter by p-value <0.05 & FC>2
      Batch_filter(
        Batch_infofile = combination,
        contentfile = list,
        outfile_suffixes = '_p005fc2',
        filepath="./",
        colname2 = 'name',
        prefixes = '',
        suffixes1 = suffixes,
        suffixes2 = suffixes2,
        P_value = 0.05,
        is.FC = 1,
        FC = 2 ,
        envelope = 1,
        is.adjustP = 0,
        is.mean_sel=0,
        mean_sel=6
      )
      
      # Filter by p-value <0.05 & FC>2 mean7
      if(grepl("miRNA",chipType,ignore.case=T)){
        Batch_filter(
          Batch_infofile = combination,
          contentfile = list,
          outfile_suffixes = '_p005fc2mean3',
          filepath="./",
          colname2 = 'name',
          prefixes = '',
          suffixes1 = suffixes,
          suffixes2 = suffixes2,
          P_value = 0.05,
          is.FC = 1,
          FC = 2 ,
          envelope = 1,
          is.adjustP = 0,
          is.mean_sel=1,
          mean_sel=3
        )
      }else{
        Batch_filter(
          Batch_infofile = combination,
          contentfile = list,
          outfile_suffixes = '_p005fc2mean7',
          filepath="./",
          colname2 = 'name',
          prefixes = '',
          suffixes1 = suffixes,
          suffixes2 = suffixes2,
          P_value = 0.05,
          is.FC = 1,
          FC = 2 ,
          envelope = 1,
          is.adjustP = 0,
          is.mean_sel=1,
          mean_sel=7
        )
      }
      
      # Filter by P-value < 0.01
      Batch_filter(
        Batch_infofile = combination,
        contentfile = list,
        outfile_suffixes = '_p001',
        filepath="./",
        colname2 = 'name',
        prefixes = '',
        suffixes1 = suffixes,
        suffixes2 = suffixes2,
        P_value = 0.01,
        is.FC = 0,
        FC = 2 ,
        envelope = 1,
        is.adjustP = 0,
        is.mean_sel=0,
        mean_sel=6
      )
      
      # Filter by p-value <0.01 & FC>2
      Batch_filter(
        Batch_infofile = combination,
        contentfile = list,
        outfile_suffixes = '_p001fc2',
        filepath="./",
        colname2 = 'name',
        prefixes = '',
        suffixes1 = suffixes,
        suffixes2 = suffixes2,
        P_value = 0.01,
        is.FC = 1,
        FC = 2 ,
        envelope = 1,
        is.adjustP = 0,
        is.mean_sel=0,
        mean_sel=6
      )
      
      #     Batch_volcanoplot(
      #     			Batch_infofile = combination,
      #     			filepath="./",
      #     			outfile_suffixes = '',
      #     			P_value = 0.05,
      #     			FC = 2 
      #     		)
      Batch_volcanoplot(Batch_infofile=combination,pvalue_cut=0.05,fc_cut=2)
      scatterplot_NoFlag(combination = combination, list = list, datafile = datafile)
      Total_AP_selection(datafile=datafile,contentfile = list,suffixes2 = suffixes2,filepathin="./",filepathout="./")
      ######  hclust PCA cor plot
      file=list.files(pattern="AP_selection.xls")[1]
      print("hclustering....!")
      hclust_plot(datafile=file,contentfile = list,suffixes2=suffixes,filepathin="./",filepathout="./",picName=picName)
      print("PCA plotting...!")
      PCA_plot(datafile=file,contentfile = list,suffixes2=suffixes,filepathin="./",filepathout="./",picName=picName)
      #prettyPCAPlot(datafile=file,list=list,suffix="_NS",picName=picName)
      print("Cor plotting...!")
      cor_plot(datafile=file,contentfile = list,colname1 = colname1,suffixes2='_NS',filepathin="./",filepathout="./")
      
      ######  heatmap plot after filter
      f_names <- list.files(pattern="VS_.+_p.+txt$")
      print("heatmapping...!")
      prettyheatmap_plot(datafile=f_names,list = list)
      #f_names <- list.files(pattern="Overlap")
      
  }#End of ttestSub
 
  #If it`s a lnc Chip!
  if(grepl("lnc.*|linc.*",chipType,ignore.case=T) ){
    if(exists("dataG")){ rm(dataG,envir=.GlobalEnv)}
    if(length(Sys.glob("Book*mRNA*.txt"))!=0 && length(Sys.glob("Book*lnc*.txt"))!=0){
      dir.create("mRNA",showWarnings=F)
      dir.create("lncRNA",showWarnings=F)
      #file.rename(from="Book_mRNA.txt",to="mRNA/Book_mRNA.txt")
      #file.rename(from="Book_lncRNA.txt",to="lncRNA/Book_lncRNA.txt")
      currentDir=getwd()
      setwd(paste(currentDir,"/mRNA",sep=""))
      ttestSub("../Book_mRNA.txt")
      setwd(paste(currentDir,"/lncRNA",sep=""))
      if(exists("dataG")){ rm(dataG,envir=.GlobalEnv)}
      ttestSub("../Book_lncRNA.txt")
      setwd(currentDir)
      
    }else{
      stop("If you r analysing lncRNA, you should have both Book_mRNA.txt and Book_lnc.txt right here!") 
    }
  }else{
    print(chipType)
    cat("Class of datafile: ",class(datafile),"\n")
    ttestSub(datafile)
    
  }  
    
    }else{
  #else if chipType is hta
  Batch_ttestHTA(datafile=datafile,
                 contentfile=list,
                 Batch_infofile=combination,
                 
                 filepathin="./",
                 filepathout="./",
                 colname1 = 'Probe_Set_ID',
                 colname2='name',
                 prefixes='',
                 suffixes='_NS',
                 T.type = ttype,
                 is.log = 1,
                 cutoff = 1,
                 is.adjustP = 0
  )
  
  
  # Filter by P-value < 0.05
  Batch_filterHTA(
    Batch_infofile = combination,
    contentfile = list,
    outfile_suffixes = '_p005',
    filepath="./",
    colname2 = 'name',
    prefixes = '',
    suffixes1 = '_NS',
    
    P_value = 0.05,
    is.FC = 0,
    FC = 2 ,
    envelope = 1,
    is.adjustP = 0,
    is.mean_sel=0,
    mean_sel=6
  )
  
  
  # Filter by p-value <0.05 & FC>2
  Batch_filterHTA(
    Batch_infofile = combination,
    contentfile = list,
    outfile_suffixes = '_p005fc2',
    filepath="./",
    colname2 = 'name',
    prefixes = '',
    suffixes1 = '_NS',
    
    P_value = 0.05,
    is.FC = 1,
    FC = 2 ,
    envelope = 1,
    is.adjustP = 0,
    is.mean_sel=0,
    mean_sel=6
  )
  
  # Filter by p-value <0.05 & FC>2 mean7
  Batch_filterHTA(
    Batch_infofile = combination,
    contentfile = list,
    outfile_suffixes = '_p005fc2mean7',
    filepath="./",
    colname2 = 'name',
    prefixes = '',
    suffixes1 = '_NS',
    
    P_value = 0.05,
    is.FC = 1,
    FC = 2 ,
    envelope = 1,
    is.adjustP = 0,
    is.mean_sel=1,
    mean_sel=7
  )
  
  # Filter by P-value < 0.01
  Batch_filterHTA(
    Batch_infofile = combination,
    contentfile = list,
    outfile_suffixes = '_p001',
    filepath="./",
    colname2 = 'name',
    prefixes = '',
    suffixes1 = '_NS',
    
    P_value = 0.01,
    is.FC = 0,
    FC = 2 ,
    envelope = 1,
    is.adjustP = 0,
    is.mean_sel=0,
    mean_sel=6
  )
  
  # Filter by p-value <0.01 & FC>2
  Batch_filterHTA(
    Batch_infofile = combination,
    contentfile = list,
    outfile_suffixes = '_p001fc2',
    filepath="./",
    colname2 = 'name',
    prefixes = '',
    suffixes1 = '_NS',
    
    P_value = 0.01,
    is.FC = 1,
    FC = 2 ,
    envelope = 1,
    is.adjustP = 0,
    is.mean_sel=0,
    mean_sel=6
  )
  
  # Batch create the volcano plot
  ###############################################################################
  # Arguments:
  # Batch_infofile : the result file of the Choose_group.r. Default is "Condition.txt"
  # filepath : the path of the input files. Default is current directory ''
  # outfile_suffixes : the suffixes of outputfile.Default is NULL
  #  P_value : the cut-off of the pvalue in T test result. Default is 0.05
  #	FC : the cutoff of the foldchange.Default is 2
  #   Batch_volcanoplot(
  #     Batch_infofile = combination,
  #     filepath="./",
  #     outfile_suffixes = '',
  #     P_value = 0.05,
  #     FC = 2 
  #   )
  Batch_volcanoplot(Batch_infofile=combination,pvalue_cut=0.05,fc_cut=2)
  # Batch create the heatmap and PCA plot
  ###############################################################################
  # Arguments:
  # datafile : the datafile name to be plot
  # contentfile : group information of samples
  # filepathin :  input file path
  # filepathout : output file path
  ######  AP selection   outputfile   _AP_selection.xls
  
  ######  hclust PCA cor plot
  hclust_plot(datafile=datafile,contentfile = list,suffixes2='_NS',filepathin="./",filepathout="./",picName=picName)
  #prettyPCAPlot(datafile=datafile,list=list,suffix="_NS",picName=picName)
  PCA_plot(datafile=datafile,contentfile = list,suffixes2='_NS',filepathin="./",filepathout="./",picName=picName)
  cor_plot(datafile=datafile,contentfile = list,colname1 = 'Probe_Set_ID',suffixes2='_NS',filepathin="./",filepathout="./")
  #######  diff_up_down_expressed_gene_count  
  f_names <- list.files(pattern="VS.+txt$")
  diff_up_down_gene_count(datafile=f_names,filepathin="./",filepathout="./")                          
  
  ######  heatmap plot after filter
  f_names <- list.files(pattern="VS_.+_p.+txt$")
  print("heatmapping...!")
  prettyheatmap_plot(datafile=f_names,list = list)
  #f_names <- list.files(pattern="Overlap")
}#End of hta 

}
