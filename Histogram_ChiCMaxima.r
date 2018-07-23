#!/bin/env Rscript

options(warn=-1)
options(scipen=999)



#The input format: 

#bait_ID	bait_chr	bait_start	bait_end	bait_name	OE_ID	OE_chr	OE_start	OE_end	OE_name	 Nreads	



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

#Function of Distance distribution between Replicates
DistanceDistribution=function(gene,file1,file2) {

  InteractionsRep1=file1[file1$Bait_name==gene,]
  InteractionsRep2=file2[file2$Bait_name==gene,]
	peaks_rep1_bed=cbind(InteractionsRep1[,c("chr_OE","start_OE","end_OE")],rep("*",length(InteractionsRep1[,"start_OE"])))
	peaks_rep2_bed=cbind(InteractionsRep2[,c("chr_OE","start_OE","end_OE")],rep("*",length(InteractionsRep2[,"start_OE"])))
	colnames(peaks_rep1_bed)=c("chr","start","end","strand")
	colnames(peaks_rep2_bed)=c("chr","start","end","strand")
		if (dim(peaks_rep1_bed)[1]>=dim(peaks_rep2_bed)[1] ) {
			f=distanceToNearest(GRanges(peaks_rep1_bed),GRanges(peaks_rep2_bed))
			f=as.data.frame(f)
  	} else if (dim(peaks_rep1_bed)[1]<dim(peaks_rep2_bed)[1]) {

			f=distanceToNearest(GRanges(peaks_rep2_bed),GRanges(peaks_rep1_bed))
			f=as.data.frame(f)
				
		}
		
			return(f)
		
}

#Function wrapper for DistanceDistribution
wraperDistanceDistribution=function(i) {

	Distances=DistanceDistribution(gene=Genes[i],file1=FirstReplicate,file2=SecondReplicate)
	return(Distances)
  
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
  cat("PromoMaxima pipeline for analyzing Promoter-Capture data: 
Usage:
    ./PromoMaxima.R [options] 
  Options:
    -f/--first			    <string>	[default:input]
    -s/--second			    <string>	[default:input]
    -o/--output         <string>  [default:Histogram_Interactions_Replicates.png]
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

  optArgs=getopt(
	rbind(
	    c('first','f', 1, 'character',"Interaction_Chr1_mEScs_Rep1.FirstReplicate"),
	    c('second', 's', 1, 'character',"Interaction_Chr1_mEScs_Rep2.FirstReplicate"),
      c('output','o',1,'character',"Histogram_Interactions_Replicates.png")
	      )
	)
}


#################################################################################
# ARGUMENTS
#################################################################################

 file1=optArgs$first
 file2=optArgs$second
 outputFile=optArgs$output

###################################################################################
# MAIN
###################################################################################




#------Main function-------#


#Reading Input Files

 FirstReplicate=read.delim(file1,header=F)
 SecondReplicate=read.delim(file2,header=F)
 #Needs a line of code to verify if there is a header or NOT and add the appropriate header
 colnames(FirstReplicate)=c("ID_Bait","chr_Bait","start_Bait","end_Bait","Bait_name","ID_OE","chr_OE","start_OE","end_OE","OE_name","N")
 colnames(SecondReplicate)=c("ID_Bait","chr_Bait","start_Bait","end_Bait","Bait_name","ID_OE","chr_OE","start_OE","end_OE","OE_name","N")
 if(class(FirstReplicate$Bait_name) =="factor")  {FirstReplicate$Bait_name=as.character(levels(FirstReplicate$Bait_name))[FirstReplicate$Bait_name]}
	if(class(FirstReplicate$chr_Bait) =="factor")   {FirstReplicate$chr_Bait=as.character(levels(FirstReplicate$chr_Bait))[FirstReplicate$chr_Bait]}
	if(class(FirstReplicate$chr_OE) =="factor")     {FirstReplicate$chr_OE=as.character(levels(FirstReplicate$chr_OE))[FirstReplicate$chr_OE]}
	if(class(FirstReplicate$OE_name) =="factor")    {FirstReplicate$OE_name=as.character(levels(FirstReplicate$OE_name))[FirstReplicate$OE_name]}
	if(class(FirstReplicate$start_OE)=="factor")    {FirstReplicate$start_OE=as.numeric(levels(FirstReplicate$start_OE))[FirstReplicate$start_OE]}
  if(class(FirstReplicate$end_OE)=="factor")       {FirstReplicate$end_OE=as.numeric(levels(FirstReplicate$end_OE))[FirstReplicate$end_OE]}
  if(class(FirstReplicate$start_Bait)=="factor")    {FirstReplicate$start_Bait=as.numeric(levels(FirstReplicate$start_Bait))[FirstReplicate$start_Bait]}
  if(class(FirstReplicate$end_Bait)=="factor")       {FirstReplicate$end_Bait=as.numeric(levels(FirstReplicate$end_Bait))[FirstReplicate$end_Bait]}
if(class(SecondReplicate$Bait_name) =="factor")  {SecondReplicate$Bait_name=as.character(levels(SecondReplicate$Bait_name))[SecondReplicate$Bait_name]}
	if(class(SecondReplicate$chr_Bait) =="factor")   {SecondReplicate$chr_Bait=as.character(levels(SecondReplicate$chr_Bait))[SecondReplicate$chr_Bait]}
	if(class(SecondReplicate$chr_OE) =="factor")     {SecondReplicate$chr_OE=as.character(levels(SecondReplicate$chr_OE))[SecondReplicate$chr_OE]}
	if(class(SecondReplicate$OE_name) =="factor")    {SecondReplicate$OE_name=as.character(levels(SecondReplicate$OE_name))[SecondReplicate$OE_name]}
	if(class(SecondReplicate$start_OE)=="factor")    {SecondReplicate$start_OE=as.numeric(levels(SecondReplicate$start_OE))[SecondReplicate$start_OE]}
  if(class(SecondReplicate$end_OE)=="factor")       {SecondReplicate$end_OE=as.numeric(levels(SecondReplicate$end_OE))[SecondReplicate$end_OE]}
  if(class(SecondReplicate$start_Bait)=="factor")    {SecondReplicate$start_Bait=as.numeric(levels(SecondReplicate$start_Bait))[SecondReplicate$start_Bait]}
  if(class(SecondReplicate$end_Bait)=="factor")       {SecondReplicate$end_Bait=as.numeric(levels(SecondReplicate$end_Bait))[SecondReplicate$end_Bait]}

 GeneReplicate1=unique(FirstReplicate$Bait_name)
 GeneReplicate2=unique(SecondReplicate$Bait_name)
 if(length(GeneReplicate1) >= length(GeneReplicate2)) { 
   Genes=GeneReplicate1 
   }else {
     Genes=GeneReplicate2
   }
#Think about cases where the gene has interactions detected in one Replicate NOT the other 

  DistanceBetweenReplicates=lapply(1:length(Genes),wraperDistanceDistribution)
  png(outputFile)
  hist(as.numeric(unlist(DistanceBetweenReplicates)))
  dev.off()
  




