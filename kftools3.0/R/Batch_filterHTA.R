Batch_filterHTA <-
function(
	Batch_infofile = 'Condition.txt',
	contentfile = 'Content.txt',
	outfile_suffixes = '_p005',
	filepath= '',
	colname2 = 'filename',
	prefixes = '',
	suffixes1 = '_NS',
	
	P_value = 0.05,
	is.FC = 1,
	FC = 2 ,
	envelope = 1,
	is.adjustP = 0,
	is.mean_sel = 1,
	mean_sel = 6
	)
{
	#-------------Read data--------------#
	#read Batch information file
	    Conditions="";
        if(class(Batch_infofile)=="data.frame"){Conditions=Batch_infofile}else{
                                     Conditions <- read.delim( Batch_infofile, sep = "\t", header = T, fill = T,check.names=FALSE)
                                                                                
    } 
    
    
    Batch_info<-subset(Conditions,Selected=='Y')
	#print 
	print(paste("The number of iterations is ",as.character(nrow(Batch_info))))
		
	#read sample file
	sampleinfo=""
         if(class(contentfile)=="data.frame"){ sampleinfo=contentfile }else{
                             sampleinfo <- read.delim(contentfile, sep = "\t", header = T, fill = T,check.names=FALSE)
    } 
	for (i in 1:nrow(Batch_info)){
		
		#read T-test result file
		## outputdir<-paste(filepath,Batch_info[i,1],"_VS_",Batch_info[i,2],sep='')
	  outputdir <- filepath 
      filename<-paste(Batch_info[i,1],"_VS_",Batch_info[i,2],sep='')
      data=""
      if(exists(paste(Batch_info[i,1],"_VS_",Batch_info[i,2],sep=''))){
      	data=get(filename)
      }else{
		data<-read.table(paste(outputdir,"/",filename,'.txt',sep=''),sep="\t",quote="",header=T,fill=T)
		assign(paste(Batch_info[i,1],"_VS_",Batch_info[i,2],sep=''),data,envir=.GlobalEnv)
	  }	
		#read label file
		label1_raw <-subset(sampleinfo,Group==as.character(Batch_info[i,1]),select=colname2)
		label1_raw<-as.matrix(label1_raw)
	 #label1<-paste(prefixes,label1_raw,suffixes,sep="")
		
		label2_raw<-subset(sampleinfo,Group==as.character(Batch_info[i,2]),select=colname2)
		label2_raw<-as.matrix(label2_raw)
	 #label2<-paste(prefixes,label2_raw,suffixes,sep="")
		
		#Extract Flag file
# 		label1_Flag<-paste(prefixes,label1_raw,suffixes2,sep="");
# 		label2_Flag<-paste(prefixes,label2_raw,suffixes2,sep="");
# 		
		#print group
		print(paste("the NO is",as.character(i)))
		print(as.character(Batch_info[i,1]))
		print(label1_raw)
		print(as.character(Batch_info[i,2]))
		print(label2_raw)

		#------------Filter procedure------------#
# 		data1_Flag<-subset(data,select=label1_Flag)#data1_Flag is a dataframe with all label1_Flag
# 		data2_Flag<-subset(data,select=label2_Flag)
# 		
# 		#Stat Flag condition of two group
# 		Flag1<-stat_Flag(data1_Flag)
# 		Flag2<-stat_Flag(data2_Flag)
# 		#Flag1 and 2 are the percentage of P for probes(lines)

		#stat_Flag(Flagdata)

# stat_Flag<-function(Flagdata)
# {	
	
# 	stat_Flag=numeric(nrow(Flagdata))
	##numeric: generate a number of 0s
# 	for (i in 1:nrow(Flagdata)){
# 		sum_PM = 0;
# 		for (j in 1:ncol(Flagdata)){
# 			if (Flagdata[i,j]!='A') sum_PM = sum_PM+1;
# 		}
# 		stat_Flag[i] = sum_PM/ncol(Flagdata);
		##stat_Flag[i] is the percentage of P in line i
# 	}

# 	#-------------output stat_Flag-----------------#
# 	stat_Flag;
	
# }


		Filter<-data
		
		
# 		#Filter the Flag
# 		Filter<-Filter_Flag(result,Propotion=envelope)
# 		Filter_Flag<-function(data,Propotion=1)
# {	
	
# 	Is=numeric(nrow(data))
	
# 	for (i in 1:nrow(data)){
# 		if (data$foldchange[i]>=1 & data$Flag1[i]>=Propotion)
# 			Is[i]=1
# 		else if (data$foldchange[i]<1&data$Flag2[i]>=Propotion) 
# 			Is[i]=1
# 		}
# 	data <- cbind(data,Is)	
# 	new <- subset(data,Is==1,select=-Is)
# }


		#Filter the P-values < 0.05(default)
		if (is.adjustP == 0){#not using the BH method
		Filter <- subset(Filter,pvalues<=P_value)
		} else {
			if ("fdr" %in% colnames(result)){ 
				Filter <- subset(Filter,fdr<=P_value)
			} else {
				stop("The data don't compute djusted p-values for simple multiple testing procedures")
			}
			
		}
		
		#Filter the foldchange
		if (is.FC==1){
			Filter <- subset(Filter,(foldchange>=FC|foldchange<=1/FC))
		} 
		
		#Filter <- subset(Filter,select=-c(Flag1,Flag2))		
		#Filter the low intensity  
		if (is.mean_sel==1){ 
		flags_tb <-  Filter[,grep(suffixes1,colnames(Filter))]
    colnames(flags_tb) <- gsub(suffixes1,"",colnames(flags_tb)) 
		group_tb_total <- sampleinfo
		group_tb <- group_tb_total[ match(colnames(flags_tb),group_tb_total[,"name"]),]
	  total_group <- as.character(unique(group_tb[,"Group"]))
		all_flags <- rep(FALSE,nrow(flags_tb))
    for(group_i in total_group){
        sample_i  <-  as.character(group_tb[grep(paste("^",group_i,"$",sep=""), group_tb[,"Group"],perl=T),"name"])
        tmp_tb <- flags_tb[,match(sample_i, colnames(flags_tb))]
        ##
        tmp_flags <- apply(tmp_tb,1,function(x){mean(x)>= mean_sel  }  )
        all_flags <- all_flags | tmp_flags 
      }
		Filter <- subset(Filter,all_flags)		
			} 
		##-------------output Filter result-----------------#
		write.table(Filter,file=paste(outputdir,"/",Batch_info[i,1],"_VS_",Batch_info[i,2],outfile_suffixes,".xls",sep=''),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
		write.table(Filter,file=paste(outputdir,"/",Batch_info[i,1],"_VS_",Batch_info[i,2],outfile_suffixes,".txt",sep=''),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)	  
  }
}
