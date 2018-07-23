#!/bin/env Rscript

options(warn=-1)
options(scipen=999)



#The input format: 





#-----------Functions-----------------------------


#Function to retrieve arguments and options from the R script 

RetreiveAndDestroy=function(opt,root,stem,regexp,SearchNames,Out,isValue,NbFound,StockNames) {
  Bool=lapply(paste(root,SearchNames,stem,sep=""),grepl,opt)
  names(Bool)=StockNames
  Pos=lapply(Bool,which)
  names(Pos)=StockNames
  disable=c()
  for(i in StockNames) {
    nbmatch=length(Pos[[i]])
    if(nbmatch>0) {
      NbFound[[i]]=NbFound[[i]]+nbmatch
      disable=c(disable,-1*Pos[[i]])
      if(is.null(Out[[i]])) {
	if(isValue[[i]]!=0) {
	  if(regexp=="next") {
	    Out[[i]]=opt[Pos[[i]]+1]
	    disable=c(disable,-1*(Pos[[i]]+1))
	  }
	  else {
	    Out[[i]]=sub(regexp,"\\1",opt[Pos[[i]]])
	  }
	}
	else {
	  Out[[i]]=TRUE
	}
      }
      else {
	if(isValue[[i]]!=0) {
	  if(regexp=="next") {
	    Out[[i]]=c(Out[[i]],opt[Pos[[i]]+1])
	    disable=c(disable,-1*(Pos[[i]]+1))
	  }
	  else {
	    Out[[i]]=c(Out[[i]],sub(regexp,"\\1",opt[Pos[[i]]]))
	  }
	}
	else {
	  Out[[i]]=c(Out[[i]],TRUE)
	}
      }
    }
  }
  if(length(disable)>0) {
    opt=opt[disable]
  }
  Out[["ARGUMENT"]]=list()
  Out[["ARGUMENT"]][["opt"]]=opt
  Out[["ARGUMENT"]][["NbFound"]]=NbFound
  return(Out)
}

getopt=function (spec=NULL,opt=commandArgs()) {
  FindArgs=which(opt=="--args")
  if(length(FindArgs)!=1) {
    stop(length(FindArgs)," --args found where 1 expected.",call.=F)
  }
  ExecName=sub("--file=","",opt[FindArgs-1])
  
  if(FindArgs<length(opt)) {
    opt=opt[(FindArgs+1):length(opt)]
  }
  else {
    opt=""
  }
  
  min.columns=5
  colNames=c("LongName","ShortName","Flag","Mod","Default")
  max.columns=6
  DimSpec=dim(spec)
  if(DimSpec[2]>min.columns) {
    colNames=c(colNames,"Description")
  }
  
  if(is.null(spec) | !is.matrix(spec) | (DimSpec[2]<min.columns | DimSpec[2]>max.columns)) {
    stop('argument "spec" is required and must be a matrix with 4|5 columns.',call.=F)
  }
  colnames(spec)=colNames
  
  spec=as.data.frame(spec,stringsAsFactors=F)
  #spec validation
  if(length(unique(c(spec$ShortName,"ARGUMENT","args")))!=DimSpec[1]+2 | length(unique(spec$LongName))!=DimSpec[1]) {
    stop('Long|Short names for flags must be unique (Long name : "ARGUMENT" and "args" forbidden).',
	"\n","List of duplicated :",
	"\n","Short: ",paste(spec$ShortName[duplicated(c(spec$ShortName,"ARGUMENT","args"))],collapse=" "),
	"\n","Long:  ",paste(spec$ShortName[duplicated(spec$LongName)],collapse=" "),call.=F)
  }
  if(length(which(nchar(spec$ShortName)>1))!=0) {
    stop('Short names flags can\'t be longer than 1 character.')
  }
  
  #initialize 
  Out=list()
  Short2Long=list()
  NbFound=list()
  isValue=list()
  for(i in 1:DimSpec[1]) {
    Short2Long[[spec$ShortName[i]]]=spec$LongName[i]
    NbFound[[spec$LongName[i]]]=0
    isValue[[spec$LongName[i]]]=spec$Flag[i]
  }
  
  #Map, retreive and suppress ARGUMENTs and arguments
  #Value ARGUMENT --example=value
  Out=RetreiveAndDestroy(opt,"^--","=.+$",".+=(.+)$",spec$LongName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #boolean ARGUMENT --example
  Out=RetreiveAndDestroy(opt,"^--","$","$",spec$LongName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #short name ARGUMENT -t value OR boolean -t
  Out=RetreiveAndDestroy(opt,"^-","$","next",spec$ShortName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #Warn about non mapped ARGUMENTs
  if(length(opt)>0) {
    PosUnkArg=which(grepl("^-",opt))
    if(length(PosUnkArg)) {
      message("Error, argument unreconized :","\n",paste(opt[PosUnkArg],collapse="\n"),"\n\n")
    }
    if(length(PosUnkArg)>0) {
      opt=opt[PosUnkArg*-1]
    }
  }
  #Arguments
  Out[["ARGUMENT"]]=opt
  
  #Validation of ARGUMENTs
  for(i in 1:DimSpec[1]) {
    if(spec$Flag[i]=="0") {#verify boolean arguments
      NbValue=length(Out[[spec$LongName[i]]])
      if(NbValue>1) {
	message("Warning : ",spec$LongName[i]," found ",NbValue," times")
      }
    }
    if(length(Out[[spec$LongName[i]]])==0) {
      Out[[spec$LongName[i]]]=spec$Default[i]
    }
    library("methods")
    Out[[spec$LongName[i]]]=as(Out[[spec$LongName[i]]],spec$Mod[i])
  }
  
  return(Out)
}



#Geometric_mean function to identify background level for each bait

Geometric_mean=function(table=table) {
	start=seq(0,1480000,by=20000)
	end=seq(20000,1500000,by=20000)
	distance=GRanges(seqnames=Rle(table$chr_Bait[1]),ranges=IRanges(start,end),strand=rep("*",length(start)))
	distance_d=as.data.frame(distance)
	distance_d$RefBinMean=1
	
	for(j in 1:length(start(distance))){
			
		bins=table[abs(table$end_OE-table$start_Bait[1]) <= end(distance)[j] & abs(table$end_OE-table$start_Bait[1]) > start(distance)[j] ,]
		distance_d$RefBinMean[j]=geometric.mean(bins[,"N"])
		if (is.na(distance_d$RefBinMean[j]) & length(distance_d$RefBinMean[j-1])!=0 ) {distance_d$RefBinMean[j]=distance_d$RefBinMean[j-1]}
		distance_d$midpoint[j]=round((end(distance)[j]+start(distance)[j])/2)
	}
  if(length(distance_d[distance_d$RefBinMean==0,]$RefBinMean)<=length(distance_d[distance_d$RefBinMean!=0,]$RefBinMean)){
	  function_distance <- glm.nb( distance_d$RefBinMean ~ distance_d$start,link="log") #Negative Binomial regression 
	  distance_d$predicted=fitted(function_distance)
   } else {
     distance_d$predicted=distance_d$RefBinMean
   }
	return (distance_d)

}

#Function local maxima To find peaks
argmax <- function(table=table,w=50,span=0.05,d=0) {
  y=table[,"N"]
	x=table[,"start_OE"]
	n <- length(y)
	y.smooth <- loess(y ~ x, span=span)$fitted
	y.max <- rollapply(zoo(y.smooth), 2*w+1,max,align="center")
	delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
	i.max <- which(delta <= d) + w
	i=0
	final=data.frame()
	distance_d=Geometric_mean(table)
	for (j in 1:length(i.max)){
		peak=table[i.max[j],]
    
		#peak=peak[!duplicated(peak$start_OE),]	
		dist=abs(peak$start_OE-peak$start_Bait[1])
		value=distance_d[dist >distance_d$start & dist <= distance_d$end,]
		value=value[!is.na(value$predicted),]
		if ( dim(value)[1] !=0 & dim(peak)[1] !=0 ) {
			if ( peak[,"N"]  > value$predicted ) {
				final[j,1]=peak$start_OE
				final[j,2]=value$predicted
				#final=final[!is.na(final[,1]),]	
			}
			
		}
		
	}
 if( dim(final)[1] !=0) {
   Final_coordinates=table[table$start_OE%in%final[,1] & table$N!=0,]
   return(Final_coordinates)
  }
}

#Main Function of identifying interactions 
IdentifyPeaks=function(gene=gene,table=ibed,window_size=50,loess_span=0.05) {	
	
  table=table[table$Bait_name==gene,]
	table=table[order(table$start_OE),]
	ref=table$start_Bait[1]
	table=table[abs(table$end_OE -ref) <= 1500000,]
	cat(paste(gene,"\n"))
  if (length(table$start_OE) !=0 & length(table$start_OE) >= 41) {
	
    	peaks=argmax(table=table,w=window_size,span=loess_span,d=0)
	    return (peaks)
 
	}
}

#Function wrapper for IdentifyPeaks
wrapeIdentifyPeaks=function(i) {

	Peaks_coordinates=IdentifyPeaks(gene=genes[i],table=ibed,window_size=50,loess_span=0.05)
	cat (paste(i,"\n"))
	return(Peaks_coordinates)
  
}




#
### =========================================================================
### ChiCMaxima's Main version beta 0.9
### -------------------------------------------------------------------------
###

if(version$major !="3" & version$minor !="0.2") {
  warning("You need 3.0.2 R Version to run the program")
}
arg=commandArgs(TRUE)
if(length(arg)==0) {
  cat("ChiCMaxima pipeline for analyzing Promoter-Capture data: 
Usage:
    ./PromoMaxima.R [options] 
  Options:
    -i/--input			    <string>	[default:input]
    -o/--output                     <string>    [default:./ouput]
    -w/--window_size                <string>    [default: 50    ]
    -s/--loess_span                 <string>    [default: 0.05  ]
   \n\n")
  q("no")

}else {
  tmp=suppressPackageStartupMessages(require("Rsamtools"))
  if(!tmp) {
    stop("Package Rsamtools required !!")
  }
  tmp=suppressPackageStartupMessages(require("GenomicRanges"))
  if(!tmp) {
    stop("Package GenomicRanges required !!")
  }
  tmp=suppressPackageStartupMessages(require("limma"))
  if (!tmp) {
	stop("Package limma required !!")
  }
  tmp=suppressPackageStartupMessages(require("MASS"))
  if (!tmp) {
	stop("Package MASS required !!")
  }
  tmp=suppressPackageStartupMessages(require("caTools"))
  if (!tmp) {
	stop("Package caTools required !!")
  }
  tmp=suppressPackageStartupMessages(require("data.table"))	
  if (!tmp) {
	stop("Package data.table required !!")
  }
  tmp=suppressPackageStartupMessages(require("base"))
  if (!tmp) {
	stop("Package base required !!")
  }
  tmp=suppressPackageStartupMessages(require("zoo"))
  if (!tmp) {
	stop("Package zoo required !!")
  }
  tmp=suppressPackageStartupMessages(require("RcppRoll"))
  if (!tmp) {
	stop("Package RcppRoll required !!")
  }
  tmp=suppressPackageStartupMessages(require("psych"))
  if (!tmp) {
	stop("Package psych required !!")
  }
  tmp=suppressPackageStartupMessages(require("plyr"))
  if(!tmp) {
	stop("Package plyr required !!")
  }
  tmp=suppressPackageStartupMessages(require("ICSNP"))
  if(!tmp) {
	stop("Package ICSNP required !!")
  }
  tmp=suppressPackageStartupMessages(require("marray")) #We can substitute this package should write our function of output
  if(!tmp){
  stop("Package marray required !!")
  }
 
  optArgs=getopt(
	rbind(
	    c('output','o', 1, 'character',"output.ibed"),
	    c('window_size', 'w', 1, 'integer', 50),
	    c('loess_span','s',1,'integer',0.05),
	    c('input','i',1,'character',"chr1_mESCs.ibed")
	    )
	)
}


#################################################################################
# ARGUMENTS
#################################################################################
files=optArgs$input
output.files=optArgs$output
window_size=optArgs$`window_size`
loess_span=optArgs$`loess_span`





###################################################################################
# MAIN
###################################################################################




#------Main function-------#

#Function to read input files 
readInput=function(file) {

	table=as.data.frame(fread(file,header=F,fill=TRUE))
 	
 	return(table)
}

doit=function(i) {

   try(readInput(files[[i]]))
  }
  data=lapply(1:length(files),doit)
  ctrl=unlist(lapply(data,class))=="try-error"

  if(sum(ctrl)>0) {
    for(i in which(ctrl)) {
     write("Problem with file :",files[i],"\n")
    }
    stop("At least a file has problem")
  }



#Convert any factor column to the right class 

for ( i in 1:length(data)){
	ibed=data[[i]]
  colnames(ibed)=c("ID_Bait","chr_Bait","start_Bait","end_Bait","Bait_name","ID_OE","chr_OE","start_OE","end_OE","OE_name","N") #Needs an if to detect a header or No Header
	if(class(ibed$Bait_name) =="factor")  {ibed$Bait_name=as.character(levels(ibed$Bait_name))[ibed$Bait_name]}
	if(class(ibed$chr_Bait) =="factor")   {ibed$chr_Bait=as.character(levels(ibed$chr_Bait))[ibed$chr_Bait]}
	if(class(ibed$chr_OE) =="factor")     {ibed$chr_OE=as.character(levels(ibed$chr_OE))[ibed$chr_OE]}
	if(class(ibed$OE_name) =="factor")    {ibed$OE_name=as.character(levels(ibed$OE_name))[ibed$OE_name]}
	if(class(ibed$start_OE)=="factor")    {ibed$start_OE=as.numeric(levels(ibed$start_OE))[ibed$start_OE]}
  if(class(ibed$end_OE)=="factor")       {ibed$end_OE=as.numeric(levels(ibed$end_OE))[ibed$end_OE]}
  if(class(ibed$start_Bait)=="factor")    {ibed$start_Bait=as.numeric(levels(ibed$start_Bait))[ibed$start_Bait]}
  if(class(ibed$end_Bait)=="factor")       {ibed$end_Bait=as.numeric(levels(ibed$end_Bait))[ibed$end_Bait]}

	genes=unique(ibed$Bait_name)
	genes=genes[order(genes)]
		
	#Identifying Peaks
	Peaks_final=lapply(1:length(genes),wrapeIdentifyPeaks)
	#ctrl=unlist(lapply(Peaks_final,class))=="try-error  #This should work! 
	#if(sum(ctrl)>0) {
    	
      		#cat("Problem with genes\n")
    			
    		#stop("At least one gene has a problem")
  	#}
  
  write.list(Peaks_final, filename = output.files[[i]])
	
} 
				




