#***********************************************************************
#*                                                                     *
#*  R code written for the R package "lmomRFA"                         *
#*                                                                     *
#*  J. R. M. HOSKING                                                   *
#*  IBM RESEARCH DIVISION                                              *
#*  T. J. WATSON RESEARCH CENTER                                       *
#*  YORKTOWN HEIGHTS                                                   *
#*  NEW YORK 10598, U.S.A.                                             *
#*                                                                     *
#*  Version 0.8    November 2008                                       *
#*                                                                     *
#***********************************************************************

cluagg<-function(x,method="ward") {
  dx<-dist(x)
  if (pmatch(method,"ward",nomatch=0)==1) dx<-dx^2
  hc<-hclust(dx,method=method)

  hcm<-(-hc$merge)
  for (j in 1:nrow(hcm)) {
    if (hcm[j,1]<0) hcm[j,1]<-hcm[-hcm[j,1],1]
    if (hcm[j,2]<0) hcm[j,2]<-hcm[-hcm[j,2],1]
    hcm[j,]<-sort(hcm[j,])
  }
  return(list(merge=hcm,wgss=cumsum(hc$height/2)))
}

cluinf<-function(merge,nclust) {
  if (is.list(merge) && names(merge)==c("merge","wgss")) merge<-merge$merge
  n<-nrow(merge)+1
  vec<-1:n
  for (i in 1:(n-nclust)) vec[vec==merge[i,2] ] <- vec[merge[i,1] ]
  assign<-match(vec,sort(unique(vec)))
  list<-lapply(1:nclust,function(i) which(assign==i))
  padnum<-formatC(1:nclust,width=nchar(nclust),format="d",flag="0")
  names(list)<-paste("cluster.",padnum,sep="")
  num<-sapply(list,length)
  return(list(assign=assign,list=list,num=num))
}

clukm<-function(x,assign,maxit=10,algorithm="Hartigan-Wong") {
  x<-as.matrix(x)
  if (nrow(x)!=length(assign))
    stop("number of rows of 'x' and length of 'assign' must be equal")
  centers<-apply(x,2,function(y) tapply(y,assign,mean))
  stats::kmeans(x,centers,iter.max=maxit,algorithm=algorithm)
}

reglmr<-function(xmom, weight) {
## Regional weighted average of L-moments
  xmom<-as.matrix(xmom)
  if (missing(weight)) weight<-rep(1,nrow(xmom))
  if (length(weight)!=nrow(xmom))
    stop("number of rows of 'xmom' and length of 'weight' must be equal")
  if (ncol(xmom)>1) xmom[,2]<-xmom[,2]/xmom[,1]
  xmom[,1]<-1
  apply(xmom,2,weighted.mean,w=weight,na.rm=TRUE)
}

regsamlmu<-function(x, nmom=5, sort.data=TRUE, lcv=TRUE) {
  if (is.list(x)) {
    xmom<-t(sapply(x,samlmu,nmom=nmom,sort.data=sort.data))
    n<-sapply(x,function(y) length(y[!is.na(y)]))
    name<-names(x)
    if (is.null(name)) name<-seq_along(x)
  } else {
    x<-as.matrix(x)
    xmom<-t(apply(x,2,samlmu,nmom=nmom,sort.data=sort.data))
    n<-apply(x,2,function(y) length(y[!is.na(y)]))
    name<-colnames(x)
    if (is.null(name)) name<-seq_len(ncol(x))
  }
  if (nmom==0) xmom<-matrix(nrow=length(name),ncol=0)
  if (nmom==1) xmom<-matrix(xmom,ncol=1,dimnames=list(NULL,"l_1"))
  if (lcv && nmom>=2) {
    xmom[,2]<-xmom[,2]/xmom[,1]
    colnames(xmom)[2]<-"t"
  }
  df<-cbind(data.frame(name=name),n,xmom)
  row.names(df)<-NULL
  return(df)
}

regtst<-function(regdat, nsim=1000){
##  Discordancy, heterogeneity and goodness-of-fit statistics
##  for regional frequency analysis

  if (!(is.data.frame(regdat) && ncol(regdat)>=7))
    stop("'regdat' must be a data frame with at least 7 columns")
  sitenames<-regdat[,1]
  if (any(duplicated(sitenames[!is.na(sitenames)])))
    warning("site names are not all different")
  if (any(regdat[,2]<=0))
    stop("record lengths must all be greater than zero")
  if (any(regdat[,3]<=0))
    stop("site means must all be greater than zero")
  if (any(regdat[,4]<0 | regdat[,4]>1))
    stop("L-CV values are not all between 0 and 1")
  if (any(regdat[,5:7]<(-1) | regdat[,5:7]>1))
    stop("L-moment ratios are not all between -1 and +1")

  nsites<-nrow(regdat)
  len<-regdat[,2]
  xmom<-t(regdat[,3:7])
  maxrec<-max(len)

  fort<-.Fortran("regtst",PACKAGE="lmomRFA",
    nsites=as.integer(nsites),
    len=as.integer(len),
    xmom=as.double(xmom),
    nsim=as.integer(nsim),
    rmom=double(5),
    d=double(nsites),
    vobs=double(3),
    vbar=double(3),
    vsd=double(3),
    h=double(3),
    z=double(5),
    para=double(30),
    rpara=double(4),
    t4fit=double(5),
    work=double(nsites*3),
    x=double(maxrec),
    maxrec=as.integer(maxrec))

  if (all(fort$d==0)) {
    is.na(d[])<-TRUE
    warning("unable to invert sum-of-squares matrix - D statistics not calculated")
  }

  if (nsim<=1) fort[c("rpara","vobs","vbar","vsd","h","z","t4fit")]<-NULL

  dc1<-c(3,3,3,3,1.3330,1.6481,1.9166,2.1401,2.3287,2.4906,
    2.6321,2.7573,2.8694,2.9709,3,3,3,3)
  dc2<-c(4,4,4,4,1.3333,1.6648,1.9821,2.2728,2.5337,2.7666,
    2.9748,3.1620,3.3310,3.4844,3.6246,3.7532,3.8718,3.9816)
  Dcrit <- if (nsites>length(dc1)) c(3,4) else c(dc1[nsites],dc2[nsites])

  para=list(
    glo=fort$para[1:3],
    gev=fort$para[6:8],
    gno=fort$para[11:13],
    pe3=fort$para[16:18],
    gpa=fort$para[21:23],
    wak=fort$para[26:30])
  names(para$glo)<-lmom:::lmom.dist$glo$parnames
  names(para$gev)<-lmom:::lmom.dist$gev$parnames
  names(para$gno)<-lmom:::lmom.dist$gno$parnames
  names(para$pe3)<-lmom:::lmom.dist$pe3$parnames
  names(para$gpa)<-lmom:::lmom.dist$gpa$parnames
  names(para$wak)<-lmom:::lmom.dist$wak$parnames

  out<-list(
    data=regdat[1:7],
    nsim=nsim,
    D=fort$d,
    Dcrit=Dcrit,
    rmom=fort$rmom,
    rpara=fort$rpara,
    vobs=fort$vobs,
    vbar=fort$vbar,
    vsd=fort$vsd,
    H=fort$h,
    para=para,
    t4fit=fort$t4fit,
    Z=fort$z)

  names(out$data)<-c("name","n","mean","t","t_3","t_4","t_5")
  names(out$rmom)<-c("mean","t","t_3","t_4","t_5")
  if (nsim>1) {
    names(out$rpara)<-lmom:::lmom.dist$kap$parnames
    names(out$t4fit)<-names(out$Z)<-c("glo","gev","gno","pe3","gpa")
  }

  class(out)<-"regtst"
  return(out)

}

print.regtst<-function(x,...) {
## Print method for an object of class "regtst"
  cat("Discordancy measures (critical value ",formatC(x$Dcrit[1],2,format="f"),")\n",sep="")
  cat(formatC(x$D,2,format="f"),"\n\n")
  if (x$nsim<=1) {
    cat("Heterogeneity measures not calculated\n")
    cat("Goodness-of-fit measures not calculated\n")
  } else {
    cat("Heterogeneity measures (based on",x$nsim,"simulations)\n")
    cat(formatC(x$H,2,format="f"),"\n\n")
    cat("Goodness-of-fit measures (based on",x$nsim,"simulations)\n")
    print(round(x$Z,2))
  }
  return(invisible(x))
}

summary.regtst<-function(object,
  prob=c(0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99,0.999),
  conf=0.90, decimals=c(4,4,2,3), ...){
## Summary method for an object of class "regtst"

  prob<-prob[prob>=0 & prob<=1]

  if (is.null(object$para)) {
    out<-c(object,list(prob=prob,quant=NULL,decimals=decimals))
  } else {
    quant<-matrix(nrow=6,ncol=length(prob))
    if (ncol(quant)>0) {
      quant[1,]<-quaglo(prob,object$para$glo)
      quant[2,]<-quagev(prob,object$para$gev)
      quant[3,]<-quagno(prob,object$para$gno)
      quant[4,]<-quape3(prob,object$para$pe3)
      quant[5,]<-quagpa(prob,object$para$gpa)
      quant[6,]<-quawak(prob,object$para$wak)
      colnames(quant)<-format(prob,scientific=FALSE)
    }
    rownames(quant)<-names(object$para)
    out<-c(object,list(conf=conf,prob=prob,quant=quant,decimals=decimals))
  }

  class(out)<-"summary.regtst"
  return(out)

}

print.summary.regtst<-function(x, decimals, ...) {
## Print an object of class "summary.regtst"
  if (missing(decimals)) decimals<-x$decimals
  ldec<-length(decimals)
  if (ldec<4) decimals<-c(decimals,c(4,4,2,3)[(ldec+1):4])
  dlmom<-decimals[1]
  dpara<-decimals[2]
  dtest<-decimals[3]
  dquant<-decimals[4]

  dat<-x$data
  dat[,3:7]<-format(dat[,3:7],digits=1,nsmall=dlmom,scientific=FALSE)
  nsites<-length(x$D)
  stars<-rep("  ",nsites)
  substring(stars,1,1)[x$D>=x$Dcrit[1] ]<-"*"
  substring(stars,2,2)[x$D>=x$Dcrit[2] ]<-"*"
  print(cbind(dat,"D(i)"=format(x$D,digits=1,nsmall=dtest),"  "=stars),
    right=FALSE)

  cat("\n")
  rmommat<-matrix(x$rmom[-1],nrow=1,
    dimnames=list("Weighted means  ",names(x$rmom)[-1]))
  print(noquote(formatC(rmommat,digits=dlmom,format="f")))

  if (any(stars!="  "))
    cat("\nFlagged test values:",
      formatC(sort(x$D[stars!="  "],decreasing=TRUE),digits=2,format="f"))

  if (x$nsim>1) {
    cat("\n")
    parmat<-matrix(x$rpara,nrow=1,
      dimnames=list("Parameters of regional kappa distribution  ",names(x$rpara)))
    print(noquote(formatC(parmat,digits=dpara,format="f")))
  }

  dnames<-c(
    "Gen. logistic      ",
    "Gen. extreme value ",
    "Gen. normal        ",
    "Pearson type III   ",
    "Gen. Pareto        ",
    "Wakeby             ")

  conf<-NA

  if (x$nsim>1) {

    if (nsites>1) {
      vobs<-formatC(x$vobs,digits=dlmom,width=dlmom+3,format="f")
      vbar<-formatC(x$vbar,digits=dlmom,width=dlmom+3,format="f")
      vsd <-formatC(x$vsd ,digits=dlmom,width=dlmom+3,format="f")
      H   <-formatC(x$H   ,digits=dtest,width=dtest+3,format="f")
      stars<-paste(format("",width=max(0,dlmom-dtest)),
        ifelse(x$H>2,"**",ifelse(x$H>1,"* ","  ")),sep="")
      cat("\n\n*****  HETEROGENEITY MEASURES  *****\n")
      cat("Number of simulations =",x$nsim,"\n\n")
      cat("Observed     s.d. of group L-CV             =",vobs[1],"\n")
      cat("Sim. mean of s.d. of group L-CV             =",vbar[1],"\n")
      cat("Sim. s.d. of s.d. of group L-CV             =",vsd[1] ,"\n")
      cat("Heterogeneity measure H[1]                  =",H[1],stars[1],"\n\n")
      cat("Observed     s.d. of L-CV / L-skew distance =",vobs[2],"\n")
      cat("Sim. mean of s.d. of L-CV / L-skew distance =",vbar[2],"\n")
      cat("Sim. s.d. of s.d. of L-CV / L-skew distance =",vsd[2] ,"\n")
      cat("Heterogeneity measure H[2]                  =",H[2],stars[2],"\n\n")
      cat("Observed     s.d. of L-skew/L-kurt distance =",vobs[3],"\n")
      cat("Sim. mean of s.d. of L-skew/L-kurt distance =",vbar[3],"\n")
      cat("Sim. s.d. of s.d. of L-skew/L-kurt distance =",vsd[3] ,"\n")
      cat("Heterogeneity measure H[3]                  =",H[3],stars[3],"\n\n")
    }

    stars<-rep(" ",5)
    if (is.numeric(x$conf) && length(x$conf)==1 && x$conf>0 && x$conf<1) {
      conf<-format(x$conf,nsmall=2)
      Zcrit<-qnorm(0.5+x$conf/2)
      stars[abs(x$Z)<=Zcrit]<-"*"
    }

    t4fit<-formatC(x$t4fit,digits=dlmom,width=dlmom+3,format="f")
    Z    <-formatC(x$Z    ,digits=dtest,width=dtest+4,format="f")
    cat("\n*****  GOODNESS-OF-FIT MEASURES  *****\n")
    cat("Number of simulations =",x$nsim,"\n\n")
    for (j in 1:5)
      cat(dnames[j],"  L-kurtosis =",t4fit[j],"   Z value =",Z[j],stars[j],"\n")
    cat("\n")

  }

  if (is.na(conf)) {
    cat("\nPARAMETER ESTIMATES\n\n")
    ok<-1:6
  } else {
    cat("\nPARAMETER ESTIMATES FOR DISTRIBUTIONS ACCEPTED AT THE",
      conf,"LEVEL\n\n")
    ok<-c(which(stars=="*"),6)
  }

  for (j in ok) { cat(dnames[j],
    formatC(x$para[[j]],digits=dpara,width=dpara+4,format="f"),"\n")
  }

  if (ncol(x$quant)>0) {
    cat("\nQUANTILE ESTIMATES\n")
    quant<-x$quant
    rownames(quant)<-dnames
    colnames(quant)<-paste(" ",colnames(quant))
    quant<-formatC(quant,digits=dquant,width=dquant+4,format="f")
    print(noquote(quant[ok,,drop=FALSE]))
  }

  return(invisible(x))
}

