# a and b, plots in rows, species in columns, except the 1. column with the plot number
# docor TRUE/FALSE: TRUE if correlations a vs b should be computed, FALSE if among
#       species correlations are given as parameter a (b is then not used)
# nonulls is the number of species permutations
# log1 if log+1 transformation of values in a and b should be done (only for docor==TRUE)
# mincount is minimal non-zero plot count for species to be species (only for docor==TRUE)
# cmethod correlation method paste to cor function, eg. Pearson,
# coval is value of correlation that split colors in the species x species plot, eg. 0.1
# onm title for boxplot
# boxonly: TRUE if only boxplots of r and r2 should be plotted, or FALSE if among species 
#       pairs are ploted too
#OUTPUTS is : data frame with the simulation number (0 are observed data),
#            if r or r2 is used, correlations statistics, t.tests for
#            diagonal ==0, and diagonal vs extradiagonal both with t value, 
#            p values and the mean estimate
conspecor<-function(a,b, nonulls=10, mincount=1,log1=TRUE, docor=TRUE,
                    cmethod ="pearson",coval=0.1,onm="",boxonly=FALSE) {
  # get the names of variables
  cargs<-as.list(match.call())
  
  if (docor) { # correlations a vs b will be computed
    print(paste("Correlation of a =",cargs$a,"and b =",cargs$b))
    
    # apply mincount
    a<-a[,c(TRUE, colSums(a[,-c(1)]>0,na.rm=T ) >= mincount)]
    b<-b[,c(TRUE, colSums(b[,-c(1)]>0,na.rm=T ) >= mincount)]
    
    # the reduction according corresponding names
    a<-a[,names(a) %in% names(b)]
    b<-b[,names(b) %in% names(a)]
    
    # get a and b stats
    asu<-colSums(a[,-c(1)],na.rm=T)
    bsu<-colSums(b[,-c(1)],na.rm=T)
    asua<-colSums(a[,-c(1)]>0,na.rm=T)
    bsua<-colSums(b[,-c(1)]>0,na.rm=T)
    
    # sort according abundance of b
    bco<-order(bsua,decreasing = T)
    a<- a[,c(1,bco+1)]
    b<- b[,c(1,bco+1)]
    
    # the same order of names and plot numbers
    if (! (all(names(a)==names(b)) &  all(names(a$plot)==names(b$plot)) ) )  
      print( "Species names or plot number are not in the same order")
    
    # set correlation matrix to empty
    cormat <- matrix(NA,nrow = ncol(a)-1,ncol = ncol(b)-1,
                     dimnames = list(names(a)[-c(1)],names(b)[-c(1)]) )
    
    # compute correlations
    if ( log1 ) { # log+1 tranform values
      al<-log(1+ a)
      bl<-log(1+ b)
      for(x in 2:ncol(a)) for(y in 2:ncol(b)) {
        cormat[x-1,y-1] <- cor( x = al[,x], y = bl[,y], method = cmethod)
      }
    } else { # non transform values
      for(x in 2:ncol(a)) for(y in 2:ncol(b)) {
        cormat[x-1,y-1] <- cor( x = a[,x], y = b[,y], method = cmethod)
      }
    }
  } else { # correlations are given as input a
    cormat <- data.frame(a)
    row.names(cormat) <- cormat[,1]
    cormat <- cormat[,-c(1)]
    # transform to numeric
    for(i in 1:ncol(cormat)) cormat[,i] <- as.numeric(cormat[,i])
    cormat<-as.matrix(cormat)
    cargs[['a']]<-"species"
    cargs[['b']]<-"species"
  }
  cormat <- round(cormat,6)
  
  if ( any(is.na(cormat)) ) print ("In the correlation matrix there are some NAs.")
  
  cormat0 <- cormat # matrix with r
  # add transform cormat matrix by ^2
  # matrix with r2
  cormat01<-apply(cormat, 2, FUN=function(x) return( ifelse (is.na(x),NA, x*x) ))
  
  # output data frame
  consperes<-data.frame(r=NA,null=NA, plusdiag=NA,negdiag=NA, plusnondiag=NA,negnondiag=NA, 
                        zerot=NA,zerop=NA, zeromean=NA, 
                        nondiat=NA,nondiap=NA, nondiamean=NA, 
                        nonnodiamean=NA)
  rdiag<-data.frame(species=NA,val=NA)
  
  # run null models and observed data analysis
  i=0; for (i in 0:nonulls) {
    
    # always start from observed data
    cormat<-cormat0
    cormat1<-cormat01
    
    #i==0 are observed data
    if (i>0) {   
      tosa<-sample(1:nrow(cormat))
      cormat<- cormat[tosa,]
      cormat1<- cormat1[tosa,]
    }
    
    # diag vs 0
    (t0<-t.test(x=diag(cormat)))
    (t0r<-t.test(x=diag(cormat1)))
    
    # diag vs non-diag
    (t1<-t.test( x=diag(cormat), y= cormat[ col(cormat) != row(cormat)]))
    (t1r<-t.test( x=diag(cormat1), y= cormat1[ col(cormat1) != row(cormat1)]))
    
    #table(diag(cormat)<0)
    negdi<-sum(diag(cormat)<0,na.rm=T)
    negdi<-c(negdi,nrow(cormat)-negdi)
    
    consperes <- rbind(consperes,
                       c("r",i, 
                         negdi, table(cormat[ col(cormat) != row(cormat)] <0),
                         t0$statistic,t0$p.value<0.05,t0$estimate, 
                         t1$statistic, t1$p.value<0.05,t1$estimate))
    consperes <- rbind(consperes, 
                       c("r2",i, 
                         negdi, table(cormat[ col(cormat) != row(cormat)] <0),
                         t0r$statistic,t0r$p.value<0.05,t0r$estimate, 
                         t1r$statistic, t1r$p.value<0.05,t1r$estimate))
    
    # only observed data are plotted
    if(i ==0) {
      if ( docor )
        write.csv(rbind(cormat,adulsBA=asu,adultsNonZero=asua,seedl=bsu,seedlNonZero=bsua),
                  paste0("cormat",".csv"))
      
      if ( !boxonly) {
        plot(c(1:ncol(cormat))~1,type='n',xlab=cargs$a,ylab=cargs$b,main="",asp=1,
             axes=F)
        t0p<-ifelse (t0$p.value <0.01,"< 0.01",paste0("= ",round(t0$p.value,2)))
        t1p<-ifelse (t1$p.value <0.01,"< 0.01",paste0("= ",round(t1$p.value,2)))
        t0rp<-ifelse (t0r$p.value <0.01,"< 0.01",paste0("= ",round(t0r$p.value,2)))
        t1rp<-ifelse (t1r$p.value <0.01,"< 0.01",paste0("= ",round(t1r$p.value,2)))
        mtext(paste0("n=",ncol(cormat)),adj=1,side=1,line=2)
        
        for(x in 1:ncol(cormat)) {
          #indivs<-cut(cormat0[x,],breaks = seq(-0.6,0.6,by=0.05))
          points(x=rep(x,ncol(cormat)),y=1:ncol(cormat),pch=15, cex=1.8,
                 col=ifelse(cormat[x,]>0.1,"darkred", 
                            ifelse(cormat[x,]>0,"coral", 
                                   ifelse(cormat[x,]< (-0.1),"darkblue","lightblue" )))
          )
          #col=colfunc(length(unique(indivs)))[as.numeric(indivs)] )
        }
        text(x=0,y=1:nrow(cormat),labels=rownames(cormat),cex=0.6,adj=1,xpd=T,
             col="darkgreen")
        text(x=1:nrow(cormat)+0.5,y=0,labels=rownames(cormat),cex=0.6,adj=1,xpd=T,
             col="darkgreen")
        #diag
        points(x=1:nrow(cormat),y=rep(-1,nrow(cormat)),pch=15, cex=1.8,xpd=T,
               col=ifelse(diag(cormat)>coval,"darkred", 
                          ifelse(diag(cormat)>0,"coral", 
                                 ifelse(diag(cormat)< (-coval),"darkblue","lightblue" )))
        )
        text(x=-1,y=-1,"Diag.:",xpd=T)
        
        mtext(paste0("t.test r: diag VS 0 (p ",t0p,")"),adj=0,line=2)
        mtext(paste0("t.test r: diag VS non-diag (p ",t1p,")"),line=1,adj=0)
        mtext(paste0("t.test r2: diag VS 0 (p ",t0rp,")"),adj=0,line=0)
        mtext(paste0("t.test r2: diag VS non-diag (p ",t1rp,")"),line=-1,adj=0)
      }
      # violin plot
      fvio <-data.frame(val=c( diag(cormat),rep(NA,(nrow(cormat)*(nrow(cormat)-1))) ),
                        class=c(rep("Diagonal",nrow(cormat)), rep("Extradiagonal",(nrow(cormat)*(nrow(cormat)-1))  )) )
      
      fvio[fvio$class == "Extradiagonal",]$val<-cormat[ col(cormat) != row(cormat)]
      boxplot(val~class,data=fvio,main=onm, ylab="r",xlab="",cex.axis=1.3, cex.lab=1.5)
      points(x=jitter(rep(1,nrow(fvio[fvio$class =="Diagonal",])),2),y=fvio[fvio$class =="Diagonal",]$val,col="green")
      points(x=jitter(rep(2,nrow(fvio[fvio$class !="Diagonal",])),2),y=fvio[fvio$class !="Diagonal",]$val,col="blue")
      ame<-aggregate(val~class,data=fvio, FUN=mean,na.rm=T)
      points(x=c(1,2),ame$val,col="red",cex=2,pch=8)
      #text(x=c(1,2),y=rep(ame$val+0.05,2),labels=round(ame$val,3),col="red",cex=1)
      abline(h=0,col="gray",lty=2)
      
      fvio[fvio$class =="Diagonal",]
      write.csv(data.frame(species=colnames(cormat),year=fvio[fvio$class =="Diagonal",]),paste0(onm,"year.csv"),row.names=F)
      
      # violin plot
      fvio <-data.frame(val=c( diag(cormat1),rep(NA,(nrow(cormat1)*(nrow(cormat1)-1))) ),
                        class=c(rep("Diagonal",nrow(cormat1)), rep("Extradiagonal",(nrow(cormat1)*(nrow(cormat1)-1))  )) )
      
      fvio[fvio$class == "Extradiagonal",]$val<-cormat1[ col(cormat1) != row(cormat1)]
      boxplot(val~class,data=fvio,main=ifelse(boxonly,"",onm), ylab="r2",xlab="",cex.axis=1.3, cex.lab=1.5)
      points(x=jitter(rep(1,nrow(fvio[fvio$class =="Diagonal",])),2),y=fvio[fvio$class =="Diagonal",]$val,col="green")
      points(x=jitter(rep(2,nrow(fvio[fvio$class !="Diagonal",])),2),y=fvio[fvio$class !="Diagonal",]$val,col="blue")
      ame<-aggregate(val~class,data=fvio, FUN=mean,na.rm=T)
      points(x=c(1,2),ame$val,col="red",cex=2,pch=8)
      #text(x=c(1,2),y=rep(ame$val+0.05,2),labels=round(ame$val,3),col="red",cex=1)
      abline(h=0,col="gray",lty=2)
      
      #library(ggplot2)
      #p<-ggplot(fvio, aes(x=class, y=val) ) + geom_violin() + geom_boxplot()
      #p
    }
  }
  consperes<- consperes[-c(1),]
  for ( i in 2:ncol(consperes)) consperes[,i]<-as.numeric(consperes[,i])
  rqua<-quantile( consperes[ (consperes$r =="r") & (consperes$null>0), ]$nondiat,probs =c(0.025,0.975),na.rm=T)
  rval<-consperes[ (consperes$r =="r") & (consperes$null==0), ]$nondiat
  r2qua<-quantile(consperes[ (consperes$r =="r2") & (consperes$null>0), ]$nondiat,probs =c(0.025,0.975),na.rm=T)
  r2val<-consperes[ (consperes$r =="r2") & (consperes$null==0), ]$nondiat
  
  print (paste0("Observed r value ",round(rval,4)," is",
                ifelse (rval>rqua[2] | rval <rqua[1]," not "," "),
                "within 95% null model range:",paste0(round(rqua,4),collapse=" ") ))
  print (paste0("Observed r2 value ",round(r2val,4)," is",
                ifelse (r2val>r2qua[2] | r2val <r2qua[1]," not "," "),
                "within 95% null model range:",paste0(round(r2qua,4),collapse=" ") ))
  return(consperes)
}
