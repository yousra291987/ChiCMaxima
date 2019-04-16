#!/bin/env Rscript

##############################################################################
# CHi-C util to analyze replicates
# Copyright (C) (2019) Yousra Ben Zouari & Tom Sexton
#
# This program is free software; you can distribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
###############################################################################

options(warn=-1)
options(scipen=999)


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

#Function of distance distribution between replicates
DistanceDistribution=function(file1,file2,gene) {

  peaks_rep1=file1[file1$Bait_name==gene,c("chr_OE","start_OE","end_OE")]
  peaks_rep2=file2[file2$Bait_name==gene,c("chr_OE","start_OE","end_OE")]
  peaks_rep1$strand=rep("*",dim(peaks_rep1)[1])
  peaks_rep2$strand=rep("*",dim(peaks_rep2)[1])
  colnames(peaks_rep1)=c("chr","start","end","strand")
  colnames(peaks_rep2)=c("chr","start","end","strand")
		if (dim(peaks_rep1)[1]>=dim(peaks_rep2)[1] ) {
			f=distanceToNearest(GRanges(peaks_rep2),GRanges(peaks_rep1),ignore.strand=T)
			f=as.data.frame(f)
  	} else  {
			f=distanceToNearest(GRanges(peaks_rep1),GRanges(peaks_rep2),ignore.strand=T)
			f=as.data.frame(f)
				
		}
			d=f$distance
			return(d)
		
}


#
### =========================================================================
### ChiCMaxima's Replicate Analysis beta 0.9
### -------------------------------------------------------------------------
###

if(version$major !="3" & version$minor !="0.2") {
  warning("You need 3.0.2 R Version to run the program")
}
arg=commandArgs(TRUE)
if(length(arg)==0) {
  cat("ChiCMaxima pipeline for analyzing Capture-HiC data: 
Usage:
    ./ChiCMaxima_RepAnalysis.R [options] 
  Options:
    -a/--file1			    <string>	[default:input1.ibed]
    -b/--file2			    <string>	[default:input2.ibed]
    -o/--output         <string>  [default:output]
       \n\n")
  q("no")

}else {
  
  tmp=suppressPackageStartupMessages(require("GenomicRanges"))
  if(!tmp) {
    stop("Package GenomicRanges required !!")
  }

  optArgs=getopt(
	rbind(
	    c('file1','a', 1, 'character',"Interaction_Chr1_mEScs_Rep1.ibed"),
	    c('file2', 'b', 1, 'character',"Interaction_Chr1_mEScs_Rep2.ibed"),
      c('output','o',1,'character',"replicate_distances")
	      )
	)
}


#################################################################################
# ARGUMENTS
#################################################################################

 file1=optArgs$file1
 file2=optArgs$file2
 output=optArgs$output

###################################################################################
# MAIN
###################################################################################


#Reading Input Files

 FirstReplicate=read.delim(file1,header=T,stringsAsFactors=F)
 SecondReplicate=read.delim(file2,header=T,stringsAsFactors=F)
 
 GeneReplicate1=unique(FirstReplicate$Bait_name)
 GeneReplicate2=unique(SecondReplicate$Bait_name)
 Genes=intersect(GeneReplicate1,GeneReplicate2)

 distances=NULL
 for (i in 1:length(Genes)) {
	distances=c(distances,DistanceDistribution(FirstReplicate,SecondReplicate,Genes[i]))
	if(i%%100==0) {
		cat(paste(i,"genes processed...\n"))
	}
 }

  hist.fn=paste0(output,"_hist.png")
  cum.fn=paste0(output,"_ecdf.png")
  quant.fn=paste0(output,"_percentiles.txt")

  png(hist.fn)
  hist(distances)
  dev.off()
  png(cum.fn)
  plot.ecdf(distances)
  dev.off()
  quant=as.data.frame(quantile(distances,seq(0,1,0.05)))
  colnames(quant)="distance"
  write.table(quant,quant.fn,col.names=T,row.names=T,quote=F,sep="\t")
  
