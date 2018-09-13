#!/bin/env Rscript

options(warn=-1)
options(scipen=999)

#-----------Functions-----------------------------


#Function to retrieve arguments and options from the R script 

RetrieveAndDestroy=function(opt,root,stem,regexp,SearchNames,Out,isValue,NbFound,StockNames) {
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
  
  #Map, retrieve and suppress ARGUMENTs and arguments
  #Value ARGUMENT --example=value
  Out=RetrieveAndDestroy(opt,"^--","=.+$",".+=(.+)$",spec$LongName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #boolean ARGUMENT --example
  Out=RetrieveAndDestroy(opt,"^--","$","$",spec$LongName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #short name ARGUMENT -t value OR boolean -t
  Out=RetrieveAndDestroy(opt,"^-","$","next",spec$ShortName,Out,isValue,NbFound,spec$LongName)
  opt=Out[["ARGUMENT"]][["opt"]]
  NbFound=Out[["ARGUMENT"]][["NbFound"]]
  Out[["ARGUMENT"]]=NULL
  #Warn about non mapped ARGUMENTs
  if(length(opt)>0) {
    PosUnkArg=which(grepl("^-",opt))
    if(length(PosUnkArg)) {
      message("Error, argument unrecognized :","\n",paste(opt[PosUnkArg],collapse="\n"),"\n\n")
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


#
### =========================================================================
### ChiCMaxima Caller version beta 0.9
### -------------------------------------------------------------------------
###

if(version$major !="3" & version$minor !="0.2") {
  warning("Optimised for R version >= 3.0.2")
}
arg=commandArgs(TRUE)
if(length(arg)==0) {
  cat("CHiCMaxima pipeline for analyzing Capture-HiC data: 
Usage:
    ./ChiCMaxima_Collate.r [options] 
  Options:
    -k/--key			    <string>	[default:collate_key.txt]
    -o/--output                     <string>    [default:./combined.ibed]
   \n\n")
  q("no")

}else {
  
  tmp=suppressPackageStartupMessages(require("data.table"))	
  if (!tmp) {
	stop("Package data.table required !!")
  }
  
  optArgs=getopt(
	rbind(
	    c('key','k', 1, 'character',"collate_key.txt"),
	    c('output', 'o', 1, 'character', "combined.ibed")
	    )
	)
}


#################################################################################
# ARGUMENTS
#################################################################################
key=optArgs$key
output=optArgs$output


###################################################################################
# MAIN
###################################################################################

#Collate parameters
key.input=read.delim(key,header=FALSE,stringsAsFactors=FALSE)
if (dim(key.input)[2]!=2) {
	stop("Bad key format - expect [IBED]  [COLUMN_NAME]\n")
}

coll=key.input[,1]
names(coll)=paste0("N.",key.input[,2])

l=vector("list")
#Check and read all ibeds
for (i in 1:length(coll)) {
	cat(paste("Reading in",coll[i],"...\n"))
	ibed=fread(coll[i],header=TRUE,fill=TRUE)
	if(!is.numeric(as.data.frame(ibed)[1,1])) {
		ibed=fread(coll[i],header=FALSE,fill=TRUE,skip=1)
	}
	if(dim(ibed)[2]!=11) {
		stop("Wrong ibed format - ID_Bait, chr_Bait, start_Bait, end_Bait, Bait_name, ID_OE, chr_OE, start_OE, OE_name, N\n")
	}
	cols = c("ID_Bait","chr_Bait","start_Bait","end_Bait","Bait_name","ID_OE","chr_OE","start_OE","end_OE","OE_name",names(coll)[i])
	for (n in 1:length(cols)) {
		setnames(ibed,n,cols[n])
	}
	l[[i]]=ibed
}


cat(paste("Merging into one ibed -",output,"...\n"))
#Merge, coercing all NAs into 0
attr=vector("list")
for (i in 1:length(l)) {
	attr[[i]]=copy(l[[i]])
	set(attr[[i]],NULL,names(coll)[i],NULL)
	remCols=names(l[[i]])[!names(l[[i]]) %in% c("ID_Bait","ID_OE",names(coll)[i])]
	for (rc in remCols) {
		set(l[[i]],NULL,rc,NULL)
	}
	setkey(l[[i]],ID_Bait,ID_OE)
}

attr=rbindlist(attr)
setkey(attr,ID_Bait,ID_OE)
attr=unique(attr)
lmerge=Reduce(function(...) merge(...,by=c("ID_Bait","ID_OE"),all=TRUE),l)
setkey(lmerge,ID_Bait,ID_OE)
lmerge=attr[lmerge]
for (n in 1:length(coll)) {
	iNcol=names(coll)[n]
	set(lmerge,which(is.na(lmerge[[iNcol]])),iNcol,0)
}
setkey(lmerge,Bait_name,start_OE)

cat("Writing table...\n")
write.table(lmerge,output,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
