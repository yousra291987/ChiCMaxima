#!/bin/env Rscript

options(warn=-1)
options(scipen=999)

all.args=commandArgs(FALSE)
fn.arg="--file="
script.name=sub(fn.arg,"",all.args[grep(fn.arg,all.args)])

args=commandArgs(T)
if (length(args)==0) {
	cat(sprintf("usage: %s <repwindow> <max oe window> <output> <replicate1 interactions> <replicate2 interactions> [replicateX interactions...]\n",script.name))
	q("no")
}

repwindow=as.numeric(args[1])
maxoe=as.numeric(args[2])
ofile=args[3]
ifiles=args[4:length(args)]
if(length(ifiles)<3) {
	stop("Not enough replicates for this script\n")
}

tmp=suppressPackageStartupMessages(require("GenomicRanges"))
if(!tmp) {
    stop("Package GenomicRanges required !!")
}



#Merge replicates with flexible window
flexwindow = function(gene,ints,repwindow) {
	sub=list()
	for (i in 1:length(ints)) {
		tab=ints[[i]]
		tab=tab[tab$Bait_name==gene,]
		sub[[i]]=tab
	}

	#define reference (query) replicate and other (subject) replicates
	n=unlist(lapply(sub,function(x) {dim(x)[1]}))
	ref=which(n==min(n))[1]
	query=sub[[ref]]
	v=1:length(sub)
	v=v[-ref]
	subjects=sub[v]

	query.bed=cbind(query[,7:9],rep("*",dim(query)[1]))
	colnames(query.bed)=c("chr","start","end","strand")
	lookup=data.frame()
	for (j in 1:length(subjects)) {
		subject=subjects[[j]]
		subject.bed=cbind(subject[,7:9],rep("*",dim(subject)[1]))
		colnames(subject.bed)=c("chr","start","end","strand")
		d=distanceToNearest(GRanges(query.bed),GRanges(subject.bed),ignore.strand=TRUE)
		d=as.data.frame(d)
		d=d[d$distance<=repwindow,]
		if(dim(d)[1]==0) {
			return(NULL)
		}
		if(j==1) {
			lookup=d[,1:2]
			colnames(lookup)=c(ref,v[j])
		} else {
			inc=remodel[d$queryHits,1]
			lookup=lookup[lookup[,1] %in% inc,]
			if(dim(lookup)[1]==0) {
				return(NULL)
			}
			lookup=cbind(lookup,d[,2])
			colnames(lookup)[1+j]=v[j]
		}
		query.bed=query.bed[d[,1],]
		subject.bed=subject.bed[d[,2],]
		query.bed$start=pmin(query.bed$start,subject.bed$start)
		query.bed$end=pmax(query.bed$end,subject.bed$end)
		remodel=cbind(lookup[,1],query.bed)
		colnames(remodel)[1]="ID"
	}
	peak=sub[[as.numeric(colnames(lookup)[1])]][lookup[,1],c(1:5,7)]
	peak$start_OE=query.bed$start
	peak$end_OE=query.bed$end
	lookup=lookup[,order(as.numeric(colnames(lookup)))]
	Ns=as.matrix(lookup)
	for (i in seq_along(as.numeric(colnames(lookup)))) {
		N.set=sub[[i]][lookup[,i],"N"]
		Ns[,i]=N.set
		peak=cbind(peak,N.set)
		colnames(peak)[dim(peak)[2]]=paste0("N.",i)
	}
	peak$N.mean=rowMeans(Ns)
	peak=peak[peak$end_OE-peak$start_OE <= maxoe,]
	if(dim(peak)[1]==0) {
		return(NULL)
	} else {
		return(peak)
	}
}



###################################################################################
# MAIN
###################################################################################

ints=list()
for (i in 1:length(ifiles)) {
	tab=read.table(ifiles[i],header=TRUE)
	if(dim(tab)[2]!=11) {
		cat(paste("Wrong ibed format:",ifiles[i],"- expecting ID_Bait, chr_Bait, start_Bait, end_Bait, Bait_name, ID_OE, chr_OE, start_OE, end_OE, OE_name, N\n"))
	q("no")
	}
	ints[[i]]=tab
}

genes=as.character(unique(ints[[1]]$Bait_name))
for (i in 2:length(ints)) {
	genes=intersect(genes,as.character(unique(ints[[i]]$Bait_name)))
}

if(length(genes)==0) {
	stop("No baits in common between replicates!\n")
}

cat(paste("Processing",length(genes),"genes...\n"))
counter=0

final.peaks=data.frame()
for (i in 1:length(genes)) {
	peaks=flexwindow(gene=genes[i],ints=ints,repwindow=repwindow)
	if(!is.null(peaks)) {
		final.peaks=rbind(final.peaks,peaks)
	}
	counter=counter+1
	if (counter %% 100 == 0) {
		cat(paste0(counter,"\n"))
	}
}

write.table(final.peaks,ofile,col.names=T,row.names=F,quote=F,sep="\t")
cat(paste(dim(final.peaks)[1],"merged interactions kept\n"))


