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

#Merge replicates with flexible window
flexwindow = function(gene,one,two,repwindow) {
	sub.one=one[one$Bait_name==gene,]
	sub.two=two[two$Bait_name==gene,]
	if (dim(sub.one)[1] <= dim(sub.two)[1]) {
		query=sub.one
		subject=sub.two
		q=1
	} else {
		query=sub.two
		subject=sub.one
		q=2
	}
	query.bed=cbind(query[,7:9],rep("*",dim(query)[1]))
	colnames(query.bed)=c("chr","start","end","strand")
	subject.bed=cbind(subject[,7:9],rep("*",dim(subject)[1]))
	colnames(subject.bed)=c("chr","start","end","strand")

	f=distanceToNearest(GRanges(query.bed),GRanges(subject.bed),ignore.strand=TRUE)
	f=as.data.frame(f)
	f=f[f$distance <= repwindow,]

	if(dim(f)[1]==0) {
		return(NULL)
	}

	peak=query[f$queryHits,c(1:5,7:10)]
	cmp=subject[f$subjectHits,c(1:5,7:10)]
	peak$start_OE=pmin(peak$start_OE,cmp$start_OE)
	peak$end_OE=pmax(peak$end_OE,cmp$end_OE)
	if(q==1) {
		peak$N.1=query[f$queryHits,"N"]
		peak$Enrichment.1=query[f$queryHits,"Enrichment"]
		peak$N.2=subject[f$subjectHits,"N"]
		peak$Enrichment.2=subject[f$subjectHits,"Enrichment"]
	} else {
		peak$N.1=subject[f$subjectHits,"N"]
		peak$Enrichment.1=subject[f$subjectHits,"Enrichment"]
		peak$N.2=query[f$queryHits,"N"]
		peak$Enrichment.2=query[f$queryHits,"Enrichment"]
	}

	return(peak)
}



#
### =========================================================================
### ChiCMaxima's Main version beta 0.9
### -------------------------------------------------------------------------
###

if(version$major !="3" & version$minor !="0.2") {
  warning("Optimised for R version >= 3.0.2")
}
arg=commandArgs(TRUE)
if(length(arg)==0) {
  cat("CHiCMaxima pipeline for analyzing Capture-HiC data: 
Usage:
    ./ChiCMaxima_MergeRep2.r [options] 
  Options:
    -a/--onepeak		    <string>	[default:replicate1_interactions.ibed]
    -b/--twopeak                    <string>    [default:replicate2_interactions.ibed]
    -d/--repdist                  <string>    [default: 0    ]
    -o/--output	                    <string>    [default: merged_interactions.ibed  ]
   \n\n")
  q("no")

}else {
  tmp=suppressPackageStartupMessages(require("GenomicRanges"))
  if(!tmp) {
    stop("Package GenomicRanges required !!")
  }
 
  optArgs=getopt(
	rbind(
	    c('onepeak','a', 1, 'character', "replicate1_interactions.ibed"),
	    c('twopeak', 'b', 1, 'character', "replicate2_interactions.ibed"),
	    c('repdist','d',1,'numeric', 0),
	    c('output','o',1,'character', "merged_interactions.ibed")
	    )
	)
}





#################################################################################
# ARGUMENTS
#################################################################################
one.file=optArgs$onepeak
two.file=optArgs$twopeak
repwindow=optArgs$repdist
ofile=optArgs$output





###################################################################################
# MAIN
###################################################################################

one=read.table(one.file,header=TRUE)
if(dim(one)[2]!=12) {
	stop("Wrong ibed format for replicate 1 - ID_Bait, chr_Bait, start_Bait, end_Bait, Bait_name, ID_OE, chr_OE, start_OE, end_OE, OE_name, N, Enrichment\n")
}
two=read.table(two.file,header=TRUE)
if(dim(one)[2]!=12) {
	stop("Wrong ibed format for replicate 2- ID_Bait, chr_Bait, start_Bait, end_Bait, Bait_name, ID_OE, chr_OE, start_OE, end_OE, OE_name, N, Enrichment\n")
}

genes=intersect(unique(one$Bait_name),unique(two$Bait_name))

if(length(genes)==0) {
	stop("No baits in common between replicates!\n")
}

cat(paste("Processing",length(genes),"genes...\n"))
counter=0

peaks_final=data.frame()
for (i in 1:length(genes)) {
	peaks_final=rbind(peaks_final,flexwindow(gene=genes[i],one=one,two=two,repwindow=repwindow))
	counter=counter+1
	if (counter %% 100 == 0) {
		cat(paste0(counter,"\n"))
	}
}

write.table(peaks_final,ofile,col.names=T,row.names=F,quote=F,sep="\t")
cat(paste(dim(peaks_final)[1],"merged interactions kept\n"))


