library(methods)
library(lattice)

getbinarydata<-function(...){ getpythondata("readbin",...) } 
getrawdata<-function(...) { getpythondata("readraw",...) }

getpythondata<-function(action, ...) {
	dots <- list(...)
	script <- getOption("bafregress", "bafRegress.py")
	args <- paste(ifelse(names(dots)=="","",paste("--", names(dots),"=",sep="")), paste("'",unlist(dots),"'", sep=""), sep="")
	cmd <- paste("python",script,action, paste(args, collapse=" "))
	#print(paste("running",cmd));
	dd<-read.table(pipe(cmd), as.is=T, header=T)
	return(dd)
}

testsamplecontamination <- function(baf, abgeno, maf, subset=NULL, ...) {
	stopifnot(all(length(baf)==length(maf), length(baf)==length(abgeno)))
	maf[maf>.5] <- 1-maf[maf>.5]
	amaf <- ifelse(abgeno==2, -maf, maf)
	callrate <- 1-mean(abgeno==3)
	genocat <- factor(abgeno, levels=c(0,2))
	subs <- abgeno==2 | abgeno==0
	if (!is.null(subset)) subs<-subs & subset
	if (all(table(genocat[subs])>0)) {
		fit <- lm(baf~amaf+genocat, subset=subs, ...)
	} else {
		fit <- lm(baf~amaf, subset=subs, ...)
	}
	a <- c(coefficients(summary(fit))[2,], callrate, nrow(fit$model))
	names(a) <- c("estimate","stderr","tval","pval","callrate", "Nhom")
	return(a)
}

testsamplesbinary<-function(file.prefix) {
	samples<-read.table(paste(file.prefix, ".samples.txt", sep=""), as.is=T, col.names=c("sample"))
	for (i in seq_along(samples$sample)) {
		#if(i %% 50==0) { print(paste(file.prefix, samples$sample[i])) }
		dd <- getbinarydata(file.prefix, sample=samples$sample[i])
		reg <- with(dd, testsamplecontamination(BAF, ABGENO, MAF))
		if (i==1) {
			allcoef <- matrix(0, ncol=length(reg), nrow=length(samples$sample))
			colnames(allcoef)<-names(reg)
		}
		allcoef[i,] <- reg
	}
    dd <- data.frame(sample=samples$sample, allcoef)
	return(dd)
}

testsamplesraw<-function(file, options) {
	samples<-do.call(getpythondata, c(list("listsamplesraw", file), options[grep("^col", names(options))]))
	for (i in seq_along(samples$sample)) {
		#if(i %% 50==0) { print(paste(file.prefix, samples$sample[i])) }
		dd <- do.call(getrawdata, c(list(file, sample=samples$sample[i]), options))
		reg <- with(dd, testsamplecontamination(BAF, ABGENO, MAF))
		if (i==1) {
			allcoef <- matrix(0, ncol=length(reg), nrow=length(samples$sample))
			colnames(allcoef)<-names(reg)
		}
		allcoef[i,] <- reg
	}
    dd <- data.frame(sample=samples$sample, allcoef)
	return(dd)

}

testonesamplebinary<-function(file.prefix, sample) {
	dd <- getbinarydata(file.prefix, sample=sample)
	testonesample(dd, sample)
}
testonesampleraw<-function(file, sample, options) {
	dd <- do.call(getrawdata, c(list(file, sample=sample), options))
	testonesample(dd, sample)
}

testonesample<-function(data, sample) {
	expcols <- c("BAF","ABGENO","MAF") 
	if (!all(expcols %in% names(data))) {
		stop(paste("Missing required field. Found:", paste(names(data), collapse=","), 
			"Expected:", paste(expcols, collapse=",")))
	}
	reg <- with(data, testsamplecontamination(BAF, ABGENO, MAF))
	res<-as.data.frame(t(reg))
	cbind(sample=sample, res)
}

plotonesamplebinary<-function(file.prefix, sample, options) {
	dd <- do.call(getbinarydata, c(list(file.prefix, sample=sample), options))
	plotonesample(dd, sample)
}
plotonesampleraw<-function(file, sample, options) {
	dd <- do.call(getrawdata, c(list(file, sample=sample), options))
	plotonesample(dd, sample)
}

plotonesample<-function(data, sample) {
	with(data, plotsample(BAF,ABGENO,MAF, main=sample))
}

plotsample <- function(baf, abgeno, maf, subset=NULL, ...) {
	require(lattice)
	xyplot(baf~maf, groups=factor(abgeno, levels=0:3, labels=c("AA","AB","BB","--")), auto.key=T, ...)
}

optparse<-function(x=commandArgs(trailingOnly=T)) {
	#simple option parser
	flag <- grep("^--",x)
	options <- list()
	if (length(flag)>0) {
		val <- flag+1
		options <- append(options, x[val])
		names(options) <- gsub("^--","",x[flag])
		args <- x[-c(flag, val)]
	} else {
		args <- x
	}
	return(list(options=options, args=args))
}

if (!interactive()) {
	cmd <- optparse()
	writeHeader = T
	outPrefix = NA
	if(!is.null(cmd$options$py)) {
		options("bafregress"=cmd$options$py)
		cmd$options$py<-NULL
	}
	if(!is.null(cmd$options$noheader)) {
		writeHeader = F
		cmd$options$noheader<-NULL
	}
	if(!is.null(cmd$options$outprefix)) {
		outPrefix = cmd$options$outprefix
		cmd$options$outprefix<-NULL
	}
	if (length(cmd$args)>0) {
		results <- switch(cmd$args[1],
			"testsamplesbinary" = testsamplesbinary(cmd$args[2]),
			"testsamplesraw" = testsamplesraw(cmd$args[2], cmd$options),
			"testonesamplebinary" = testonesamplebinary(cmd$args[2], cmd$args[3], cmd$options),
			"testonesampleraw" = testonesampleraw(cmd$args[2], cmd$args[3], cmd$options),
			"plotsamplebinary" = plotonesamplebinary(cmd$args[2], cmd$args[3], cmd$options),
			"plotsampleraw" = plotonesampleraw(cmd$args[2], cmd$args[3], cmd$options),
			{cat(paste("Unrecognized command:",cmd$args[1],"\n")); data.frame()}
		)
		if (is(results, "data.frame") && nrow(results)>0) {
			write.table(results, stdout(), 
				col.names=writeHeader, row.names=F, quote=F, sep="\t")
		}
		if(is(results, "trellis")) {
			if (!is.na(outPrefix)) {
				png(paste(outPrefix, "png", sep="."))
			} else {
				png()
			}
			print(results)
			dev.off()
		}
	}
}
