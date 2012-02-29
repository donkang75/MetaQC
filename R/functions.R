GetPVal <- function(dat) {
	group <- factor(colnames(dat))
	if(any(is.na(dat))) {
		p.value <- apply(dat, 1, function(x) {
					t.test(x[which(group==levels(group)[1])], x[which(group==levels(group)[2])], var.equal=T)$p.value 
				}) 
	} else {
		p.value <- t.test2(dat, group)
	}
	return(p.value)
}

#From SAM
t.test2 <- function(x,y) {
	varr <- function(x, meanx=NULL){
		n <- ncol(x)
		p <- nrow(x)
		Y <-matrix(1,nrow=n,ncol=1)
		if(is.null(meanx)){   meanx <- rowMeans(x)}
		ans<- rep(1, p)
		xdif <- x - meanx %*% t(Y)
		ans <- (xdif^2) %*% rep(1/(n - 1), n)
		ans <- drop(ans)
		return(ans)
	}
	n1 <- sum(y==levels(y)[1])
	n2 <- sum(y==levels(y)[2])
	
	m1 <- rowMeans(x[,y==levels(y)[1],drop=F])
	m2 <- rowMeans(x[,y==levels(y)[2],drop=F])
	
	df <- n1+n2-2
	sd <- sqrt( ((n2-1) * varr(x[, y==levels(y)[2]], meanx=m2) + (n1-1) * varr(x[, y==levels(y)[1]], meanx=m1) )*(1/n1+1/n2)/df )
	
	tstat <-  (m2 - m1) / sd
	pval <- 2 * pt(-abs(tstat), df)
	return(pval)
}

GetEWPval <- function(PDat) {
	pchisq(-2*rowSums(log(PDat),na.rm=TRUE), df=2*rowSums(!is.na(PDat)), lower.tail=F)
}

printLog <- function(msg, verbose) {
	if(verbose)
		print(sprintf("[%s] %s", Sys.time(), msg))
}

#From gtools
combinations <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
	if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
			0) 
		stop("bad value of n")
	if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
			0) 
		stop("bad value of r")
	if (!is.atomic(v) || length(v) < n) 
		stop("v is either non-atomic or too short")
	if ((r > n) & repeats.allowed == FALSE) 
		stop("r > n and repeats.allowed=FALSE")
	if (set) {
		v <- unique(sort(v))
		if (length(v) < n) 
			stop("too few different elements")
	}
	v0 <- vector(mode(v), 0)
	if (repeats.allowed) 
		sub <- function(n, r, v) {
			if (r == 0) 
				v0
			else if (r == 1) 
				matrix(v, n, 1)
			else if (n == 1) 
				matrix(v, 1, r)
			else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
										1, r, v[-1]))
		}
	else sub <- function(n, r, v) {
			if (r == 0) 
				v0
			else if (r == 1) 
				matrix(v, n, 1)
			else if (r == n) 
				matrix(v, 1, n)
			else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
						Recall(n - 1, r, v[-1]))
		}
	sub(n, r, v[1:n])
}

#From matrixStats
rowVars <- function (x, center = NULL, ...) {
	n <- !is.na(x)
	n <- rowSums(n)
	n[n <= 1] <- NA
	if (is.null(center)) {
		center <- rowMeans(x, ...)
	}
	x <- x - center
	x <- x * x
	x <- rowSums(x, ...)
	x <- x/(n - 1)
	x
}

GMT2List <- function(f, saveAs=NULL) {
	trim <- function(x) {
		sub(" +$", "", sub("^ +", "", x))
	}
	
	con <- file(f, "r", blocking = FALSE)
	GList <- readLines(con) # empty
	GList <- lapply(GList, function(x) unlist(strsplit(x, "\t", fixed=T)))
	GList <- foreach(x=iter(GList)) %dopar% sapply(strsplit(x, "///", fixed=T), function(xx) trim(xx[[1]]))
	GList.names <- unlist(foreach(x=iter(GList)) %dopar% x[1])
	GList <- foreach(x=iter(GList)) %dopar% x[-(1:2)]
	names(GList) <- GList.names
	close(con)
	if(!is.null(saveAs))
		save(GList, file=saveAs)
	return(GList)
}

getFileExt <- function(x) {
	sub(".*[.]", "", x)
}

getFileName <- function(x) {
	sub("(.+)[.][^.]+$", "\\1", basename(x))
}

union.rec <- function(.list, ...){
	if(length(.list)==1) return(.list[[1]])
	Recall(c(list(union(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

Download <- function(pkg, fn) {
	url <- paste("http://cloud.github.com/downloads/donkang75", pkg, fn, sep="/")
	res <- try(download.file(url, fn, mode = "wb"), silent=TRUE)
	return(res)
}

samr2 <- function(data,  resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass","One class", "Two class paired","Two class unpaired timecourse",
				"One class timecourse","Two class paired timecourse", "Pattern discovery"), s0=NULL, s0.perc=NULL, nperms=100, center.arrays=FALSE, 
		testStatistic=c("standard","wilcoxon"), time.summary.type=c("slope","signed.area"), regression.method=c("standard","ranks"), 
		return.x=FALSE, knn.neighbors=10, random.seed=NULL, xl.mode=c("regular","firsttime","next20","lasttime"), xl.time=NULL,  xl.prevfit=NULL){
	
# Modified to implement parallel execution by Don Kang at Aug 11, 2010
	
# for robustness of code, we define constants to be used later
	samr.const.quantitative.response <- "Quantitative"
	samr.const.twoclass.unpaired.response <- "Two class unpaired"
	samr.const.survival.response <- "Survival"
	samr.const.multiclass.response <- "Multiclass"
	samr.const.oneclass.response <- "One class"
	samr.const.twoclass.paired.response <- "Two class paired"
	samr.const.twoclass.unpaired.timecourse.response <- "Two class unpaired timecourse"
	samr.const.oneclass.timecourse.response <- "One class timecourse"
	samr.const.twoclass.paired.timecourse.response <- "Two class paired timecourse"
	samr.const.patterndiscovery.response <- "Pattern discovery"
	
	## individual functions for each response type
	
	ttest.func <- function(x,y,s0=0, sd=NULL){
		
		n1 <- sum(y==1)
		n2 <- sum(y==2)
		
		p <- nrow(x)
		m1 <- rowMeans(x[,y==1,drop=F])
		m2 <- rowMeans(x[,y==2,drop=F])
		
		
		if(is.null(sd)){
			sd <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m2) + (n1-1) * varr(x[, y==1], meanx=m1) )*(1/n1+1/n2)/(n1+n2-2) )
		}
		
		numer <-  m2 - m1
		
		dif.obs <- (numer)/(sd + s0)
		return(list(tt=dif.obs,numer=numer, sd=sd))
	}
	
	
	wilcoxon.func <- function(x,y,s0=0){
		
		n1 <- sum(y==1)
		n2 <- sum(y==2)
		p=nrow(x)
		
		r2= rowSums(t(apply(x,1,rank))[,y==2,drop=F])
		numer=r2- (n2/2)*(n2+1) -(n1*n2)/2
		sd=sqrt(n1*n2*(n1+n2+1)/12)
		tt= (numer)/(sd + s0)
		return(list(tt=tt,numer=numer, sd=rep(sd,p)))
	}
	
	
	
	
	onesample.ttest.func <- function(x,y,s0=0, sd=NULL){
		n <- length(y)
		x <- x*matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=TRUE)
		m <- rowMeans(x)
		if(is.null(sd)){ sd <- sqrt( varr(x, meanx=m)/n)}
		dif.obs <- m/(sd + s0)
		return(list(tt=dif.obs, numer=m,sd=sd))
		
	}
	
	patterndiscovery.func=function(x,s0=0, eigengene.number=1){
		a=mysvd(x, n.components=eigengene.number)
		v=a$v[,eigengene.number]
		
# here we try to guess the most interpretable orientation for the eigengene
		om=abs(a$u[, eigengene.number]) >quantile(abs(a$u[, eigengene.number]),.95)
		if(median(a$u[om,eigengene.number])<0){
			v=-1.0*v
		}
		
		aa=quantitative.func(x,v,s0=s0)
		eigengene=cbind(1:nrow(a$v),v)
		dimnames(eigengene)=list(NULL,c("sample number","value"))
		return(list(tt=aa$tt, numer=aa$numer, sd=aa$sd, eigengene=eigengene))
		
		
	}
	
	
	
	
	paired.ttest.func <- function(x,y,s0=0, sd=NULL){
		
		nc <- ncol(x)/2
		o <- 1:nc
		o1 <- rep(0,ncol(x)/2);o2 <- o1
		for(j in 1:nc){o1[j] <- (1:ncol(x))[y==-o[j]]}
		for(j in 1:nc){o2[j] <- (1:ncol(x))[y==o[j]]}
		d <- x[,o2,drop=F]-x[,o1,drop=F]
		su <- x[,o2,drop=F]+x[,o1,drop=F]
		if(is.matrix(d)){ 
			m <-  rowMeans(d)
		}
		if(!is.matrix(d)) {m <- mean(d)}
		if(is.null(sd)){
			if(is.matrix(d)){ sd <- sqrt(varr(d, meanx=m)/nc)}
			if(!is.matrix(d)){sd <- sqrt(var(d)/nc)}
		}
		
		
		dif.obs <- m/(sd+s0)
		return(list(tt=dif.obs, numer=m, sd=sd))
		
		
	}
	
	
	cox.func <- function(x,y,censoring.status,s0=0){
		coxscor <- function(x, y, ic, offset = rep(0., length(y))) {
			
			## computes cox scor function for rows of nx by n matrix  x
			
			## first put everything in time order
			
			n <- length(y)
			nx <- nrow(x)
			yy <- y + (ic == 0.) * (1e-05)
			otag <- order(yy)
			y <- y[otag]
			ic <- ic[otag]
			x <- x[, otag, drop = F]
			
			##compute  unique failure times, d=# of deaths at each failure time, 
			##dd= expanded version of d to length n, s=sum of covariates at each
			## failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
			
			offset <- offset[otag]
			
			a <- coxstuff(x, y, ic, offset = offset)
			nf <- a$nf
			fail.times <- a$fail.times
			s <- a$s
			d <- a$d
			dd <- a$dd
			nn <- a$nn
			nno <- a$nno
			w <- rep(0., nx)
			for(i in (1.:nf)) {
				w <- w + s[, i]
				oo<- (1.:n)[y >= fail.times[i]]
				r<-rowSums(x[, oo, drop = F] * exp(offset[oo]))
				w<- w - (d[i]/nno[i])*r 
			}
			
			return(list(scor = w, coxstuff.obj = a))
			
		}
		
		coxvar <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL){
			
			## computes information elements (var) for cox
			## x is nx by n matrix of expression  values
			
			nx <- nrow(x)
			n <- length(y)
			yy <- y + (ic == 0.) * (1e-06)
			otag <- order(yy)
			y <- y[otag]
			ic <- ic[otag]
			x <- x[, otag, drop = F]
			offset <- offset[otag]
			if(is.null(coxstuff.obj)) {
				coxstuff.obj <- coxstuff(x, y, ic, offset = offset)
			}
			
			nf <- coxstuff.obj$nf
			fail.times <- coxstuff.obj$fail.times
			s <- coxstuff.obj$s
			d <- coxstuff.obj$d
			dd <- coxstuff.obj$dd
			nn <- coxstuff.obj$nn
			nno <- coxstuff.obj$nno
			
			x2<- x^2
			oo <- (1.:n)[y >= fail.times[1] ]
			sx<-(1/nno[1])*rowSums(x[, oo] * exp(offset[oo]))
			s<-(1/nno[1])*rowSums(x2[, oo] * exp(offset[oo]))
			w <-  d[1] * (s - sx * sx)
			
			
			for(i in 2.:nf) {
				oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]
				sx<-(1/nno[i])*(nno[i-1]*sx-rowSums(x[, oo,drop=F] * exp(offset[oo])))
				s<-(1/nno[i])*(nno[i-1]*s-rowSums(x2[, oo,drop=F] * exp(offset[oo])))
				w <- w + d[i] * (s - sx * sx)
			}
			return(w)
			
		}
		
		coxstuff<- function(x, y, ic, offset = rep(0., length(y))) {
			
			fail.times <- unique(y[ic == 1.])
			nf <- length(fail.times)
			n <- length(y)
			nn <- rep(0., nf)
			nno <- rep(0., nf)
			for(i in 1.:nf) {
				nn[i] <- sum(y >= fail.times[i])
				nno[i] <- sum(exp(offset)[y >= fail.times[i]])
			}
			
			s <- matrix(0., ncol = nf, nrow = nrow(x))
			d <- rep(0., nf)
			
			##expand d out to a vector of length n
			
			for(i in 1.:nf) {
				o <- (1.:n)[(y == fail.times[i]) & (ic == 1.)]
				d[i] <- length(o)
			}
			
			oo <- match(y, fail.times)
			oo[ic==0]<-NA
			oo[is.na(oo)]<- max(oo[!is.na(oo)])+1
			s<-t(rowsum(t(x),oo))
			if(ncol(s)> nf){s<-s[,-ncol(s)]}
			dd <- rep(0., n)
			
			for(j in 1.:nf) {
				dd[(y == fail.times[j]) & (ic == 1.)] <- d[j]
			}
			
			return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))
			
		}
		
		scor <- coxscor(x,y, censoring.status)$scor
		sd <- sqrt(coxvar(x,y, censoring.status))
		tt <- scor/(sd+s0)
		return(list(tt=tt, numer=scor, sd=sd))
	}
	
	multiclass.func <- function(x,y,s0=0){
		
		##assumes y is coded 1,2...
		
		nn <- table(y)
		m <- matrix(0,nrow=nrow(x),ncol=length(nn))
		v <- m
		for(j in 1:length(nn)){
			m[,j] <- rowMeans(x[,y==j], na.rm=T)
			v[,j] <- (nn[j]-1)*varr(x[,y==j], meanx=m[,j])
		}
		mbar <- rowMeans(x)
		mm <- m-matrix(mbar,nrow=length(mbar),ncol=length(nn))
		fac <- (sum(nn)/prod(nn))
		scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=TRUE)*mm*mm,1,sum)))
		
		sd <- sqrt(rowSums(v)*(1/sum(nn-1))*sum(1/nn))
		tt <- scor/(sd+s0)
		return(list(tt=tt, numer=scor, sd=sd))
	}
	
	
	
	quantitative.func  <-
			function(x,y,s0=0){
		
# regression of x on y
		
		my=mean(y)
		yy <- y-my
		temp <- x%*%yy
		mx=rowMeans(x)
		syy= sum(yy^2)
		
		scor <- temp/syy
		b0hat <- mx-scor*my
		ym=matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=T)
		xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+ym*matrix(scor,nrow=nrow(x),ncol=ncol(x))
		sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
		sd <- sigma/sqrt(syy)
		tt <- scor/(sd+s0)
		
		return(list(tt=tt, numer=scor, sd=sd))
		
	}
	
	
	est.s0<-function(tt,sd,s0.perc=seq(0,1, by=.05)){
		
		
		## estimate s0 (exchangeability) factor for denominator.
		## returns the actual estimate s0 (not a percentile)
		
		br=unique(quantile(sd,seq(0,1,len=101)))
		nbr=length(br)
		
		a<-cut(sd,br,labels=F)
		a[is.na(a)]<-1
		cv.sd<-rep(0,length(s0.perc))
		
		for(j in 1:length(s0.perc)){
			w<-quantile(sd,s0.perc[j])
			w[j==1]<-0
			tt2<-tt*sd/(sd+w)
			tt2[tt2==Inf]=NA
			sds<-rep(0,nbr-1)
			
			for(i in 1:(nbr-1)){
				sds[i]<-mad(tt2[a==i], na.rm=TRUE)
			}
			
			cv.sd[j]<-sqrt(var(sds))/mean(sds)
		}
		
		o=(1:length(s0.perc))[cv.sd==min(cv.sd)]
		
# we don;t allow taking s0.hat to be  0th percentile when min sd is 0
		s0.hat=quantile(sd[sd!=0],s0.perc[o])
		
		return(list(s0.perc=s0.perc,cv.sd=cv.sd, s0.hat= s0.hat))
		
	}
	
	
	
	
	
	
	
	
	
	
	varr <- function(xx, meanx=NULL){
		n <- ncol(xx)
		p <- nrow(xx)
		Y <-matrix(1,nrow=n,ncol=1)
		if(is.null(meanx)){   meanx <- rowMeans(xx, na.rm=T)}
		ans<- rep(1, p)
		xdif <- xx - meanx %*% t(Y)
		ans <- (xdif^2) %*% rep(1/(n - 1), n)
		ans <- drop(ans)
		return(ans)
		
	}
	
	parse.block.labels.for.2classes=function(y){
#this only  works for 2 class case- having form jBlockn, where j=1 or 2
		n=length(y)
		y.act=rep(NA,n)
		blocky=rep(NA,n)
		for(i in 1:n){
			blocky[i]=as.numeric(substring(y[i],7,nchar(y[i])))
			y.act[i]=as.numeric(substring(y[i],1,1))
		}
		return(list(y.act=y.act,blocky=blocky))
	}
	
	parse.time.labels.and.summarize.data=function(x,y,resp.type, time.summary.type){
# parse time labels, and summarize time data for each person, via a slope or area
# does some error checking too
		
		n=length(y)
		
		
		last5char=rep(NA,n)
		last3char=rep(NA,n)
		for( i in 1:n){
			last3char[i]=substring(y[i],nchar(y[i])-2, nchar(y[i]))
			last5char[i]=substring(y[i],nchar(y[i])-4, nchar(y[i]))
		}
		
		if(sum(last3char=="End")!= sum(last5char=="Start")){
			stop("Error in format of  time course data: a Start or End tag is missing")
		}
		
		y.act=rep(NA,n)
		timey=rep(NA,n)
		person.id=rep(NA,n)
		k=1
		end.flag=FALSE
		person.id[1]=1
		if(substring(y[1],nchar(y[1])-4,nchar(y[1]))!="Start"){ 
			stop("Error in format of  time course data: first cell should have a Start tag")
		}
		
		for(i in 1:n){
			cat(i)
			j=1
			while(substring(y[i],j,j)!="T"){j=j+1}
			end.of.y= j-1
			y.act[i]=as.numeric(substring(y[i],1,end.of.y))
			timey[i]=substring(y[i],end.of.y +5 ,nchar(y[i]))
			if(nchar(timey[i])>3 & substring(timey[i],nchar(timey[i])-2,nchar(timey[i]))=="End"){
				end.flag=TRUE
				timey[i]=substring(timey[i],1,nchar(timey[i])-3)
			}
			if(nchar(timey[i])>3 & substring(timey[i],nchar(timey[i])-4,nchar(timey[i]))=="Start"){  timey[i]=substring(timey[i],1,nchar(timey[i])-5)}
			
			if( i<n & !end.flag){ person.id[i+1]=k}
			if( i<n & end.flag){k=k+1; person.id[i+1]=k}
			
			end.flag=FALSE                
		}
		timey=as.numeric(timey)
		
		
		
		# do a check that the format was correct
		tt=table(person.id, y.act)
		junk=function(x){sum(x!=0)}
		if(sum(apply(tt,1,junk)!=1)>0){
			num=(1:nrow(tt))[apply(tt,1,junk)>1]
			stop(paste("Error in format of  time course data, timecourse #",as.character(num)) )
		}
		
		npeople=length(unique(person.id))
		
		newx=matrix(NA,nrow=nrow(x),ncol=npeople)
		sd=matrix(NA,nrow=nrow(x),ncol=npeople)
		
		
		for(j in 1:npeople){
			jj=person.id==j
			tim=timey[jj]
			xc=t(scale(t(x[,jj, drop=F]),center=TRUE,scale=FALSE))
			if(time.summary.type=="slope"){
				junk=quantitative.func(xc,tim-mean(tim))
				newx[,j]=junk$numer
				sd[,j]=junk$sd
			}
			if(time.summary.type=="signed.area"){
				junk=timearea.func(x[,jj, drop=F],tim)
				newx[,j]=junk$numer
				sd[,j]=junk$sd
				
			}
		}
		
		y.unique=y.act[!duplicated(person.id)]
		return(list(y=y.unique,  x=newx, sd=sd))
	}
	
	
	check.format=function(y, resp.type, censoring.status=NULL){
		
# here i do some format checks for the input data$y
# note that checks for time course data are done in the parse function for time course;
#  we then check the output from the parser in this function
		
		
		
		
		
		
		if(resp.type==samr.const.twoclass.unpaired.response | resp.type==samr.const.twoclass.unpaired.timecourse.response){
			if(sum(y==1)+sum(y==2) !=length(y)){
				stop(paste("Error in input response data: response type ",
								resp.type, " specified; values must be 1 or 2"))
			}
		}
		
		if(resp.type==samr.const.twoclass.paired.response | resp.type==samr.const.twoclass.paired.timecourse.response ){
			if(sum(y)!=0){
				stop(paste("Error in input response data: response type ",
								resp.type, " specified; values must be -1, 1, -2, 2, etc"))
			}
			if(sum(table(y[y>0])!=abs(table(y[y<0]))) ){
				stop(paste("Error in input response data:  response type ",
								resp.type, " specified; values must be -1, 1, -2, 2, etc"))
			}
		}
		
		
		
		if(resp.type==samr.const.oneclass.response | resp.type==samr.const.oneclass.timecourse.response){
			if(sum(y==1) !=length(y)){
				stop(paste("Error in input response data: response type ",
								resp.type, " specified;  values must all be 1"))
			}
		}
		
		
		
		
		if(resp.type==samr.const.multiclass.response){
			tt=table(y)
			nc=length(tt)
			if(sum(y<=nc & y> 0) <length(y)){
				stop(paste("Error in input response data: response type ",
								resp.type, " specified; values must be 1,2, ... number of classes"))
			}
			
			for(k in 1:nc){
				if(sum(y==k)<2){
					stop(paste("Error in input response data: response type ",
									resp.type, " specified; there must be >1 sample per class"))
				}}
		}
		
		
		if(resp.type==samr.const.quantitative.response){
			if(!is.numeric(y)){
				stop(paste("Error in input response data: response type", resp.type, " specified; values must be numeric"))
			}
		}
		
		if(resp.type==samr.const.survival.response){
			if(is.null(censoring.status)){
				stop(paste("Error in input response data: response type ",
								resp.type, " specified; error in censoring indicator"))
			}
			if(!is.numeric(y) | sum(y<0) >0){
				stop(paste("Error in input response data:  response type ",resp.type, " specified; survival times  must be numeric and nonnegative"))
				if(sum(censoring.status==0) +sum(censoring.status==1) !=length(censoring.status)){
					stop(paste("Error in input response data: response type ",
									resp.type, " specified; censoring indicators must be 0 (censored) or 1 (failed)"))
					
				}
			}
			if(sum(censoring.status==1)<1){
				stop(paste("Error in input response data:   response type ",resp.type, " specified; there are no uncensored observations"))
			}
		}
		
		return()
	}
	
	mysvd<-function(x,  n.components=NULL){
# finds PCs of matrix x
		p<-nrow(x)
		n<-ncol(x)
		
# center the observations (rows)
		
		feature.means<-rowMeans(x)
		x<- t(scale(t(x),center=feature.means,scale=F))
		
		
		if(is.null(n.components)){n.components=min(n,p)}
		if(p>n){
			a<-eigen(t(x)%*%x)
			v<-a$vec[,1:n.components,drop=FALSE]
			d<-sqrt(a$val[1: n.components,drop=FALSE])
			
			u<-scale(x%*%v,center=FALSE,scale=d)
			
			
			return(list(u=u,d=d,v=v))
		}
		else{
			
			junk<-svd(x,LINPACK=TRUE)
			nc=min(ncol(junk$u), n.components)
			return(list(u=junk$u[,1:nc],d=junk$d[1:nc],
							v=junk$v[,1:nc]))
		}
	}
	
	permute.rows <-function(x)
	{
		dd <- dim(x)
		n <- dd[1]
		p <- dd[2]
		mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
		matrix(t(x)[order(mm)], n, p, byrow = TRUE)
	}
	
	this.call=match.call()  
	resp.type.arg=match.arg(resp.type)
	
	
	xl.mode=match.arg(xl.mode)
	
	if(!is.null(random.seed)){
		set.seed(random.seed)
	}
	
	
# initialize some things (harmlessly), just so that xl.mode will work correctly
	
	
	x=NULL
	censoring.status=NULL
	
	testStatistic <- match.arg(testStatistic)
	time.summary.type <- match.arg(time.summary.type)
	regression.method <- match.arg(regression.method)
	
	x=data$x
	y=data$y
	argy=y
	if(!is.null(data$eigengene.number)){
		eigengene.number=data$eigengene.number
	}
	
	# impute missing data
	if(sum(is.na(x))>0){
		require(impute)
		x=impute.knn(x,k=knn.neighbors)
		if(!is.matrix(x)){x=x$data}
	}
	
	are.blocks.specified=FALSE
	
	# center columns of  array data if requested
	if(center.arrays){
		x<-scale(x,center=apply(x,2,median),scale=FALSE)
	}
	
	# check if there are blocks for 2 class unpaired case
	if(resp.type==samr.const.twoclass.unpaired.response){
		if(substring(y[1],2,6)=="Block" | substring(y[1],2,6)=="block"){
			junk=parse.block.labels.for.2classes(y)
			y=junk$y; blocky=junk$blocky
			are.blocks.specified=TRUE
		}}
	
	# make sure 1,2, -1,1,, etc are non-character values  coming from Excel
	if( resp.type==samr.const.twoclass.unpaired.response |resp.type==samr.const.twoclass.paired.response | resp.type==samr.const.oneclass.response | resp.type==samr.const.quantitative.response | resp.type==samr.const.multiclass.response) {y=as.numeric(y)}
	# parse and summarize, if timecourse data
	
	sd.internal=NULL
	
	if(resp.type==samr.const.twoclass.unpaired.timecourse.response | 
			resp.type==samr.const.twoclass.paired.timecourse.response |
			resp.type==samr.const.oneclass.timecourse.response){
		junk=parse.time.labels.and.summarize.data(x,y, resp.type, time.summary.type)
		y=junk$y;
		x=junk$x;
		sd.internal=sqrt(rowMeans(junk$sd^2))
		if(min(table(y))==1){
			cat("",fill=T)
			cat("Warning: only one timecourse in one or more classes;
							SAM plot and FDRs will be unreliable; only gene scores are informative",fill=T)
		}
	}
	
	# if the data is timecourse, we have already summarized the time aspect.
	# Thus we change the resp.type to the appropriate non-time-course type. Note that the original value
	#  of resp.type was saved above in resp.type.arg
	
	if(resp.type==samr.const.twoclass.unpaired.timecourse.response){ resp.type=samr.const.twoclass.unpaired.response}
	if(resp.type==samr.const.twoclass.paired.timecourse.response){ resp.type=samr.const.twoclass.paired.response}
	if(resp.type==samr.const.oneclass.timecourse.response){ resp.type=samr.const.oneclass.response}
	
	stand.contrasts=NULL
	stand.contrasts.95=NULL
	
	if(resp.type==samr.const.survival.response){censoring.status=data$censoring.status}
	
	# do a thorough error  checking of the response data
	check.format(y,resp.type=resp.type,censoring.status=censoring.status)
	
	# transform to ranks if appropriate
	if(resp.type==samr.const.quantitative.response & regression.method=="ranks")
	{
		y=rank(y)  
		x=t(apply(x,1,rank))
	}
	
	n <- nrow(x)
	ny <- length(y)
	sd <- NULL;numer <- NULL
	
	# initial computation to get sd
	if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="standard"){init.fit <- ttest.func(x,y, sd=sd.internal);numer <- init.fit$numer;sd <- init.fit$sd}
	if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="wilcoxon"){init.fit <- wilcoxon.func(x,y);numer <- init.fit$numer;sd <- init.fit$sd}
	
	if(resp.type==samr.const.oneclass.response){init.fit <- onesample.ttest.func(x,y, sd=sd.internal);numer <- init.fit$numer;sd <- init.fit$sd}
	if(resp.type==samr.const.twoclass.paired.response){init.fit <- paired.ttest.func(x,y, sd=sd.internal); numer <- init.fit$numer ; sd <- init.fit$sd}
	
	if(resp.type==samr.const.survival.response){init.fit <- cox.func(x,y,censoring.status);numer <- init.fit$numer;sd <- init.fit$sd}
	if(resp.type==samr.const.multiclass.response){init.fit <- multiclass.func(x,y); numer <- init.fit$numer;sd <- init.fit$sd}
	if(resp.type==samr.const.quantitative.response){init.fit <- quantitative.func(x,y);numer <- init.fit$numer;sd <- init.fit$sd}
	if(resp.type==samr.const.patterndiscovery.response){init.fit <- patterndiscovery.func(x);numer <- init.fit$numer;sd <- init.fit$sd}
	
	# for wilcoxon or rank regression or patterndiscovery , we set s0 to the 5th percentile of the sd values (automatic
	# estimation is not possible as the values of sd are too coarse)
	
	# also if dataset is small (< 500 genes), we don't attempt automatic
	# estimation of s0
	
	if( (resp.type==samr.const.quantitative.response & (testStatistic=="wilcoxon" | regression.method=="ranks" ) | resp.type==samr.const.patterndiscovery.response) |
			resp.type==samr.const.twoclass.unpaired.response & testStatistic=="wilcoxon" | 
			(nrow(x)<500) & is.null(s0) & is.null(s0.perc))
	{s0=quantile(sd,.05); s0.perc= 0.05}
	
	# estimate s0 if necessary
	if(is.null(s0)){
		if(!is.null(s0.perc)){
			if((s0.perc != -1 & s0.perc < 0) | s0.perc > 100){
				stop("Illegal value for s0.perc: must be between 0 and 100, or equal
								to (-1) (meaning that s0 should be set to zero)")
			}
			if(s0.perc== -1){s0=0}
			if(s0.perc>=0) {s0 <- quantile(init.fit$sd,s0.perc/100)}
		}
		
		if(is.null(s0.perc)){   
			s0=est.s0(init.fit$tt,init.fit$sd)$s0.hat
			s0.perc=100*sum(init.fit$sd<s0)/length(init.fit$sd)
		}
	}
	
	
	# compute test statistics on original data
	if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="standard"){tt <- ttest.func(x,y,s0=s0, sd=sd.internal)$tt}
	if(resp.type==samr.const.twoclass.unpaired.response & testStatistic=="wilcoxon"){tt <- wilcoxon.func(x,y,s0=s0)$tt}
	
	if(resp.type==samr.const.oneclass.response ){tt <- onesample.ttest.func(x,y,s0=s0, sd=sd.internal)$tt}
	if(resp.type==samr.const.twoclass.paired.response){tt <- paired.ttest.func(x,y,s0=s0, sd=sd.internal)$tt}
	if(resp.type==samr.const.survival.response){tt <- cox.func(x,y,censoring.status,s0=s0)$tt}
	if(resp.type==samr.const.multiclass.response){
		junk2 <- multiclass.func(x,y,s0=s0)
		tt=junk2$tt
		stand.contrasts=junk2$stand.contrasts
	}
	
	if(resp.type==samr.const.quantitative.response){tt <- quantitative.func(x,y,s0=s0)$tt}
	
	if(resp.type==samr.const.patterndiscovery.response){
		junk <- patterndiscovery.func(x,s0=s0, eigengene.number=eigengene.number);
		tt<- junk$tt;
	}
	
	2*pnorm(tt, lower.tail=FALSE)
}



