options(warn=-1)
suppressPackageStartupMessages(require(tcltk2))
suppressPackageStartupMessages(require(tkrplot))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(caTools))
suppressPackageStartupMessages(require(rtracklayer))

#Note of called 'global' variables
#N - the N columns from original ibed
N=NULL
#settings - list for each condition level: name, N (names of N columns), cols (index from ibed)
settings=NULL
#ibed - the ibed table of all called reads
ibed=NULL
#set - user input for managing replicates, used to generate settings
set=list(lab=list(),entry=list())
#pset - user input for managing conditions for plotting
pset=list()
#ints - interaction sets
ints=list()
#pints - tmp store of interaction parameters
pint=list(name=list(),plot=list())
#intplots - tmp store of if interaction is plotted or not
intplots=list()
#genes track
genome=NULL
#epigenomic tracks
epigenome=list()
#trks - tmp store of track parameters
trks=list(lab=list(),entry=list())
#cols, the default colors for each plot
plotcols=c("blue","red","green","black","cyan","purple","yellow","gray")


###################################################################################
#Called interaction overlays

#Change color helper - interactions
changeintcolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=ints[[to.change]]$color, title="Choose Color"))))
	tkconfigure(manint$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		ints[[to.change]]$color <- color
		tkconfigure(manint$env$canvas[[to.change]], bg=color)
	}
}

#Interaction plot settings util
verify.int.settings=function() {
	for (i in 1:length(ints)) {
		ints[[i]]$name <- as.character(tkget(pint$name[[i]]))
		ints[[i]]$plot <- as.character(tclvalue(intplots[[i]]))
	}
	tkdestroy(manint)
}

#Interaction plot settings
manage.ints=function() {
	manint <- tktoplevel()
	tktitle(manint)="Plot settings for interactions"
	labcol=tklabel(manint,text="Color")
	labplot=tklabel(manint,text="Plot")
	manint$env$change = tk2button(manint, text="Change Color", command=changeintcolor)
	to.change <- tkentry(manint,width=5,textvariable=tclVar("1"))
	tkgrid(labcol, row=1, column=2, padx=15, pady=5)
	tkgrid(labplot, row=1, column=3, padx=15, pady=5)
	tkgrid(manint$env$change, row=1, column=4, padx=15, pady=5)
	tkgrid(to.change, row=1, column=5, padx=15, pady=5)
	for (i in 1:length(ints)) {
		pint$name[[i]] <- tkentry(manint,width=10,textvariable=tclVar(ints[[i]]$name))
		manint$env$canvas[[i]] = tk2canvas(manint,width=30,height=10,bg=ints[[i]]$color)
		manint$env$plot[[i]] <- tk2checkbutton(manint)
		intplots[[i]] <- tclVar(ints[[i]]$plot)
		tkconfigure(manint$env$plot[[i]], variable=intplots[[i]])
		tkgrid(pint$name[[i]],row=i+1,column=1,padx=15,pady=5)
		tkgrid(manint$env$canvas[[i]],row=i+1,column=2,padx=15,pady=5)
		tkgrid(manint$env$plot[[i]],row=i+1,column=3,pady=5)
	}
	done=tkbutton(manint,text="OK",command=verify.int.settings)
	tkgrid(done,row=length(ints)+2,column=1,pady=5)
}

#Load interactions
load.ints=function() {
	files=tclvalue(tkgetOpenFile(filetypes = "{{Interaction files} {.ibed}} {{All files} *}",multiple=TRUE))
	if(!nchar(files)) {
		return()
	}
	files=unlist(strsplit(files,split=" ",fixed=FALSE,perl=FALSE,useBytes=FALSE))
	for (i in 1:length(files)) {
		ints[[i]] <- list()
		int.name = unlist(strsplit(files[i],"/"))
		int.name=int.name[length(int.name)]
		int.name=unlist(strsplit(int.name,"\\."))
		int.name=paste(int.name[-(length(int.name))],collapse="")
		ints[[i]]$name <- int.name
		ints[[i]]$color <- plotcols[i]
		ints[[i]]$plot <- "1"
		ints[[i]]$ibed <- read.table(files[i],header=TRUE,stringsAsFactors=FALSE)
	}
}

##################################################################
#Genes and epigenomic tracks

#Gene track loading
load.genes=function() {
	fn = tclvalue(tkgetOpenFile(filetypes = "{{gene files} {.txt}} {{all files} *}"))
	if(!nchar(fn)) {
		return()
	}
	genome <- read.table(fn,header=TRUE,stringsAsFactors=FALSE)
}

#Epigenomic tracks

#Change color helper function - tracks
changetrkcolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=epigenome[[to.change]]$color, title="Choose Color"))))
	tkconfigure(mantrk$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		epigenome[[to.change]]$color <<- color
		tkconfigure(mantrk$env$canvas[[to.change]], bg=color)
	}
}

#Track plot setting util
verify.trk.settings=function() {
	for (i in 1:length(epigenome)) {
		level = as.integer(tkget(trks$entry[[i]]))
		epigenome[[i]]$level <- level
	}
	tkdestroy(mantrk)
}

#Plotting setting for tracks
manage.tracks=function() {
	mantrk <- tktoplevel()
	tktitle(mantrk)="Plot settings for tracks"
	labcol=tklabel(mantrk,text="Color")
	lablev=tklabel(mantrk,text="Level")
	mantrk$env$change = tk2button(mantrk,text="Change Color", command=changetrkcolor)
	to.change <- tkentry(mantrk,width=5,textvariable=tclVar("1"))
	tkgrid(labcol,row=1,column=2,padx=15,pady=5)
	tkgrid(lablev,row=1,column=3,padx=15,pady=5)
	tkgrid(mantrk$env$change,row=1,column=4,padx=15,pady=5)
	tkgrid(to.change,row=1,column=5,padx=15,pady=5)
	for (i in 1:length(epigenome)) {
		trks$lab[[i]] <- tklabel(mantrk,text=names(epigenome)[i])
		mantrk$env$canvas[[i]] = tk2canvas(mantrk,width=30,height=10,bg=epigenome[[i]]$color)
		trks$entry[[i]] <- tkentry(mantrk,width=5,textvariable=tclVar(epigenome[[i]]$level))
		tkgrid(trks$lab[[i]],row=i+1,column=1,padx=15,pady=5)
		tkgrid(mantrk$env$canvas[[i]],row=i+1,column=2,padx=15,pady=5)
		tkgrid(trks$entry[[i]],row=i+1,column=3,padx=15,pady=5)
	}
	done=tkbutton(mantrk,text="OK",command=verify.trk.settings)
	tkgrid(done,row=length(epigenome)+2,column=1,pady=5)
}

#Import bigwig/bedgraph tracks
import.bw=function(files) {
	for (i in 1:length(files)) {
		cat(paste("Reading file:",files[i],"\n"))
		trackname = unlist(strsplit(files[i],split="/",fixed=FALSE,perl=FALSE,useBytes=FALSE))
		trackname = trackname[length(trackname)]
		trackname = unlist(strsplit(trackname,split="\\."))[1]
		epigenome[[trackname]] <- list()
		epigenome[[trackname]]$level <- 0
		epigenome[[trackname]]$color <- plotcols[i]
		epigenome[[trackname]]$track <- import(files[i])
		cat("Finished reading file\n")
	}
}

#Epigenomic track load wrapper
load.epigenome=function() {
	input = tclvalue(tkgetOpenFile(filetypes = "{{Track Files} {.bw .bedGraph}} {{All files} *}",multiple=TRUE))
	if(!nchar(input)) {
		return()
	}
	files = unlist(strsplit(input,split=" ",fixed=FALSE,perl=FALSE,useBytes=FALSE))
	import.bw(files)
	manage.tracks()
}




###################################################################################
#ChiC profile inputs


#Change color helper function - condition/replicate set
changecolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=settings[[to.change]]$color, title="Choose Color"))))
	tkconfigure(manage$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		settings[[to.change]]$color <<- color
		tkconfigure(manage$env$canvas[[to.change]], bg=color)
	}
}


#Condition/replicate plot setting util
verify.plot.settings=function() {
	for (i in 1:length(settings)) {
		settings[[i]]$name <- as.character(tkget(pset[[i]]))
	}
	tkdestroy(manage)
}

#Plotting settings for condition/replicates
manage.set=function() {
	manage <- tktoplevel()
	tktitle(manage)="Plot settings for conditions"
	labname=tklabel(manage,text="Name")
	labcol=tklabel(manage,text="Color")
	manage$env$change = tk2button(manage, text="Change Color", command = changecolor)
	to.change <- tkentry(manage,width=5,textvariable=tclVar("1"))
	tkgrid(labname, row=1, column=2, padx=15, pady=5)
	tkgrid(labcol, row=1, column=3, padx=15, pady=5)
	tkgrid(manage$env$change, row=1, column=4, padx = 15, pady=5)
	tkgrid(to.change, row=1, column=5, padx=5, pady=5)
	for (i in 1:length(settings)) {
		lablevel=tklabel(manage, text=paste(settings[[i]][["N"]],collapse=","))
		tkgrid(lablevel, row=i+1, column=1, sticky="e",pady=5)
		pset[[i]] <- tkentry(manage,width=10,textvariable=tclVar(settings[[i]]$name))
		manage$env$canvas[[i]] = tk2canvas(manage, width=30, height=10, bg=settings[[i]]$color)
		tkgrid(pset[[i]],row=i+1, column=2, padx=15, pady=5)
		tkgrid(manage$env$canvas[[i]],row=i+1,column=3, padx=15, pady=5)
	}
	done=tkbutton(manage,text="OK",command=verify.plot.settings)
	tkgrid(done,row=length(settings)+2,column=1,pady=5)
}


#Condition/replicate setting util
verify.settings=function() {
	assign=NULL
	for (i in 1:length(N)) {
		assign[i] = as.integer(tkget(set$entry[[i]]))
	}
	names(assign) = N
	keep=assign[is.integer(assign) & assign>0]
	levels=unique(keep)
	settings <- list()
	for (l in 1:length(levels)) {
		settings[[levels[l]]] <- list()
		hits = which(assign==levels[l])
		names = names(hits)
		name=unlist(strsplit(names[1],"\\."))[1]
		settings[[levels[l]]][["name"]] <- name
		settings[[levels[l]]][["N"]] <- names(hits)
		settings[[levels[l]]][["cols"]] <- hits+10
		settings[[levels[l]]][["color"]] <- plotcols[l]
	}
	tkdestroy(select)
	return(settings)
}


#Identify the N columns, and allocate to condition/replicates
set.conditions=function(ibed) {
	cols=names(ibed)
	if(length(cols)<11) {
		cat("File must have first ten ibed columns, plus sequence read columns\n")
		return(NULL)
	}
	N <- cols[-(1:10)]
	if(length(N)!=length(unique(N))) {
		cat("Sequence read columns must be unique\n")
		return(NULL)
	}
	select <- tktoplevel()
	tktitle(select)="Manage conditions and replicates"
	for (i in 1:length(N)) {
		set$lab[[i]] <- tklabel(select,text=N[i])
		set$entry[[i]] <- tkentry(select,width=10,textvariable=tclVar("0"))
		tkgrid(set$lab[[i]],row=i,column=1,sticky="e",pady=5)
		tkgrid(set$entry[[i]],row=i,column=2,sticky="w",pady=5)
	}
	done=tkbutton(select,text="OK",command=verify.settings)
	tkgrid(done,row=length(N)+1,column=1,pady=5)
}

#Open ibed from dialog box, auto prompt for conditions setting
load.ibed = function() {
	ifile = tclvalue(tkgetOpenFile(filetypes = "{{ibed files} {.ibed}} {{all files} *}"))
	if(!nchar(ifile)) {
		return()
	}
	cat(paste0("opening file ",ifile,"\n"))
	ibed <- as.data.frame(fread(ifile,header=TRUE,fill=TRUE),stringsAsFactors=FALSE)
	set.conditions(ibed)
}


#############################################################################
#Save screenshot to file
save.figure=function() {
	fn=tclvalue(tkgetSaveFile(filetypes = "{{Browser screenshots} {.eps}} {{All files} *}"))
	if(nchar(fn)==0) {
		return()
	}
	cat(sprintf("saving figure to: %s\n",fn))
	setEPS()
	postscript(fn)
	gene=tclvalue(tkget(bait))
	min.plot=as.numeric(tclvalue(tkget(min.plot)))
	max.plot=as.numeric(tclvalue(tkget(max.plot)))
	win.plot=as.numeric(tclvalue(tkget(win.plot)))
	plot.chic(ibed=ibed,gene,min.plot,max.plot,win.plot,settings=settings)
	dev.off()
}

###########################################################################
#Plot profile

#Gene plot util
parking=function(left,right) {
	y=rep(-1,length(right))
	lengths=right-left
	for (i in order(lengths,decreasing=TRUE)) {
		otherleft=left
		otherleft[i]=NA
		otherright=right
		otherright[i]=NA
		placed=FALSE
		y[i]=0
		while(placed==FALSE) {
			placed=sum((right[i]>otherleft[y==y[i]])&(left[i]<otherright[y==y[i]]),na.rm=TRUE)==0
			if(placed==FALSE) {
				y[i]=y[i]+1
			}
		}
	}
	return(y)
}

#Gene plot sub
plot.genes=function(genome,x.min,x.max) {
	y_plot=parking(genome$Start,genome$End)
	plot(c(x.min,x.max),c(1,-max(y_plot)-0.5),col="white",ylab="",xlab="",fg="white",col.axis="white",xaxs="i",yaxs="i")
	arrowHeads=pretty(x.min:x.max,n=50)
	for(i in 1:dim(genome)[1]) {
		x=c(genome$Start[i],arrowHeads[arrowHeads>genome$Start[i]&arrowHeads<genome$End[i]],genome$End[i])
		if(genome$Strand[i]=="-") {
			arrows(x[2:length(x)],-y_plot[i],x[1:length(x)-1],col="blue",length=0.08)
		} else {
			arrows(x[1:length(x)-1],-y_plot[i],x[2:length(x)],col="blue",length=0.08)
		}
		text(genome$Start[i],-y_plot[i]+0.4,adj=0,labels=genome$Name[i])
	}
}

#Obtain track autoscale levels
get.track.levels=function(epigenome) {
	levels=NULL
	for (i in 1:length(epigenome)) {
		if(epigenome[[i]]$level>0) {
			levels=c(levels,epigenome[[i]]$level)
		}
	}
	levels=unique(levels[order(levels)])
	if(levels[1]!=1 & levels[length(levels)]!=length(levels)) {
		cat("Track levels must be ascending integers with no gaps: 1,2,3,...\n")
		return(NULL)
	}
	plottracklevels=list()
	for (l in 1:length(levels)) {
		for (i in 1:length(epigenome)) {
			if(epigenome[[i]]$level == levels[l]) {
				if(length(plottracklevels)<l) {
					plottracklevels[[l]]=names(epigenome)[i]
				} else {
					plottracklevels[[l]]=c(plottracklevels[[l]],names(epigenome[i]))
				}
			}
		}
	}
	return(plottracklevels)
}

#Track plot sub
plot.tracks=function(chr,x.min,x.max, tracks) {
	tmps=list()
	plotlim=numeric()
	for (i in 1:length(tracks)) {
		tmp=epigenome[[tracks[i]]]$track
		tmp=tmp[as.character(seqnames(tmp))==chr & start(tmp)>=x.min & end(tmp)<=x.max,]
		tmps[[i]]=tmp
		plotlim=max(c(plotlim,score(tmp)))
	}
	for (i in 1:length(tmps)) {
		plot.new()
		plot.window(xlim=c(x.min,x.max),ylim=c(0,plotlim),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
		segments(start(tmps[[i]]),0,end(tmps[[i]]),score(tmps[[i]]),col=epigenome[[tracks[i]]]$color)
		title(main=tracks[i],col.main=epigenome[[tracks[i]]]$color,adj=0,line=-2)
	} 
}

#Plot profile
plot.chic = function(ibed,gene,min.plot,max.plot,win.plot,settings) {
	table=ibed[ibed$Bait_name == gene,]
	table=table[order(table$start_OE),]
	ref=table$start_Bait[1]

	#Normalize across all used interaction sets
	usecols = NULL
	for (i in 1:length(settings)) {
		usecols=c(usecols,settings[[i]][["cols"]])
	}
	n=normalizeBetweenArrays(cbind(table[,usecols]))
	for (i in 1:dim(n)[2]) {
		table[,usecols[i]]=n[,i]
	}
	
	#Define plot window (max/min overrides plot.window)
	if (!is.na(min.plot) & !is.na(max.plot)) {
		table = table[table$start_OE >= min.plot & table$end_OE <= max.plot,]
	}
	else if (!is.na(win.plot)) {
		table=table[abs(table$start_OE - ref)<= as.integer(win.plot),]
	}
	else {
		cat("Must select either a window size or a full min/max constraint\n")
		return(NULL)
	}
	
	n.plots=1
	if(length(epigenome)>0) {
		for (i in 1:length(epigenome)) {
			if(epigenome[[i]]$level>0) {
				n.plots=n.plots+1
			}
		}
	}
	if(!is.null(genome)) {
		n.plots=n.plots+1
	}
	if(n.plots>1) {
		layout(matrix(1:n.plots,ncol=1,nrow=n.plots),heights=c(20,rep(10,n.plots-1)))
	}
	#first plot - CHiC
	tab=cbind(table[,settings[[1]][["cols"]]])
	avg = rowMeans(tab)
	avg = caTools::runmean(avg,k=10,endrule="mean",align="center")
	x.min = min(table$start_OE)
	x.max = max(table$end_OE)
	y.max = max(avg)
	main.text = paste0(gene,"- ",table$chr_Bait[1],": ",x.min,"-",x.max,"\n")
	plot(x=table$start_OE, y=avg,type="l",col=settings[[1]][["color"]],ylab="Running mean QN signal",xlab="Genomic coordinates",xaxs="i",yaxs="i",ylim=c(0,y.max),lwd=2,main=main.text)
	if(length(settings)>1) {
		for (i in 2:length(settings)) {
		tab=cbind(table[,settings[[i]][["cols"]]])
		avg=rowMeans(tab)
		lines(x=table$start_OE, y=caTools::runmean(avg,k=10,endrule="mean",align="center"),col=settings[[i]]$color,lwd=2)
		}
	}
	plotnames=c()
	plotcolors=c()
	for (i in 1:length(settings)) {
		plotnames=c(plotnames,settings[[i]]$name)
		plotcolors=c(plotcolors,settings[[i]]$color)
	}
	legend("topleft",inset=0.05,plotnames,col=plotcolors,lty=1.5,title="ChiC profile")
	abline(v=ref,col="black")

	#Overlay interactions
	if(length(ints)>0) {
		keep=list()
		for (i in 1:length(ints)) {
			if(ints[[i]]$plot == "1") {
				keep[[ints[[i]]$name]]=ints[[i]]
			}
		}
		if(length(keep)>0) {
			for (i in 1:length(keep)) {
				inttab=keep[[i]]$ibed
				inttab=inttab[inttab$Bait_name==gene,]
				inttab=inttab[inttab$start_OE>=x.min & inttab$end_OE<=x.max,]
				intcol=keep[[i]]$color
				rect(inttab$start_OE,rep(0,dim(inttab)[1]),inttab$end_OE,rep(y.max,dim(inttab)[1]),border=intcol)
			}
		}
		intcols=c()
		for (i in 1:length(keep)) {
			intcols=c(intcols,keep[[i]]$color)
		}
		legend("left",inset=0.05,names(keep),col=intcols,lty=1.5,title="Interactions")
	}
	n.plots=n.plots-1
	#second plot - gene track
	if(!is.null(genome)) {
		genome = genome[genome$Chr==table$chr_Bait[1] & genome$Start>=x.min & genome$End<=x.max,]
		if(dim(genome)[1]>0) {
			plot.genes(genome,x.min,x.max)
		}
	n.plots=n.plots-1
	}
	#other plots - epigenomic tracks, with level-specific autoscaling
	if(n.plots>0) {
		plottracklevels <- get.track.levels(epigenome)
		for (l in 1:length(plottracklevels)) {
			plot.tracks(table$chr_Bait,x.min,x.max,plottracklevels[[l]])
		}
	}
}

#Wrapper to set up plotting function and output to gui
refresh.gui = function() {
	if(is.null(ibed)) {
		cat("Need to open an ibed file first - see menu file\n")
		return(NULL)
	}
	gene = tclvalue(tkget(bait))
	min.plot = as.numeric(tclvalue(tkget(min.plot)))
	max.plot = as.numeric(tclvalue(tkget(max.plot)))
	win.plot = as.numeric(tclvalue(tkget(win.plot)))
	hscale = 1.99
	vscale = 1.99
	if(is.null(browse$env$plot)) {
		browse$env$plot = tkrplot(browse, plot.chic(ibed=ibed,gene,min.plot,max.plot,win.plot,settings=settings),hscale=hscale,vscale=vscale)
		tkgrid(browse$env$plot)
	}
	else {
		tkrreplot(browse$env$plot, plot.chic(ibed=ibed,gene,min.plot,max.plot,win.plot,settings=settings),hscale=hscale,vscale=vscale)
		tkgrid(browse$env$plot)
	}
}

########################################################################################################################

#GUI

main = tktoplevel()
browse = tktoplevel()
tktitle(main) = "Main"
tktitle(browse) = "Browser"
main$env$menu = tk2menu(main)
tkconfigure(main,menu=main$env$menu)

#Main menu (File -> New,Save Image,Quit)
main$env$menuFile = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="File",menu=main$env$menuFile)
tkadd(main$env$menuFile, "command", label = "New", command=load.ibed)
tkadd(main$env$menuFile, "command", label = "Save Image", command=save.figure)
tkadd(main$env$menuFile, "command", label = "Quit", command=function() tkdestroy(main))

#Conditions menu (-> Set Conditions,Plot Conditions)
main$env$menuCon = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Conditions", menu=main$env$menuCon)
tkadd(main$env$menuCon, "command", label="Set Conditions", command=function() set.conditions(ibed))
tkadd(main$env$menuCon, "command", label="Plot Conditions", command=manage.set)

#Interactions menu (-> Load Interactions,Manage Interaction Plots)
main$env$menuInt = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Interactions", menu=main$env$menuInt)
tkadd(main$env$menuInt, "command", label="Load Interactions", command=load.ints)
tkadd(main$env$menuInt, "command", label="Manage Interaction Plots", command=manage.ints)

#Tracks menu (-> Load Genes, Load Tracks, Manage Tracks)
main$env$menuTrk = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Tracks", menu=main$env$menuTrk)
tkadd(main$env$menuTrk, "command", label="Load Genes", command=load.genes)
tkadd(main$env$menuTrk, "command", label="Load Tracks", command=load.epigenome)
tkadd(main$env$menuTrk, "command", label="Manage Tracks", command=manage.tracks)


#Plot window and bait choice
labmin = tklabel(main,text="start coordinate")
min.plot <- tkentry(main, width=20, textvariable=tclVar("NA"))
tkgrid(labmin, row=1, column=1, sticky="e")
tkgrid(min.plot, row=1, column=2, sticky="w")

labmax = tklabel(main,text="end coordinate")
max.plot <- tkentry(main, width=20, textvariable=tclVar("NA"))
tkgrid(labmax, row=2, column=1, sticky="e")
tkgrid(max.plot, row=2, column=2, sticky="w")

labwin = tklabel(main,text="plot window")
win.plot <- tkentry(main, width=10, textvariable=tclVar("1000000"))
tkgrid(labwin, row=3, column=1, sticky="e")
tkgrid(win.plot, row=3, column=2, sticky="w")

labbait = tklabel(main,text="bait")
bait <- tkentry(main, width=10, textvariable=tclVar("Myc"))
tkgrid(labbait, row=4, column=1, sticky="e")
tkgrid(bait, row=4, column=2, sticky="w")

ok=tkbutton(main,text="OK",command=refresh.gui)
tkgrid(ok,row=5,column=1)
