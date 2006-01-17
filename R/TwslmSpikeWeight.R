#********************************************************************#
#                                                                    #
#   Function to calculate TW-SLM (two-way semi-linear model          #
#   for normalization of cDNA microarray data using spike data       #
#   and spots quality information                                    #
#                                                                    #
#********************************************************************#

.First.lib <- function(lib, pkg) {
    library.dynam("TwslmSpikeWeight", pkg, lib)
}

library(splines)                #load bs package

#****************************************************************
#             OLS normalization                                 *
#****************************************************************

spike<-function(sld,geneid,rt,intn,weight=NULL,s.sld,s.geneid,s.rt,s.intn,s.weight=NULL,
                s.constant=0.0,df=12,degree=3,tol=1e-5){

     cat("Calculation is running..."," ",date(),"\n")
     EPS<-tol
     id<-as.integer(unclass(factor(geneid)))     #make id from 1,2,....J
     slide<-as.integer(unclass(factor(sld)))     #make slide 1,2,...,n
     ii<-order(slide,id)
     
     id<-id[ii]
     slide<-slide[ii]
     ratio<-rt[ii]
     intensity<-intn[ii]
     uniqueid<-unique(geneid[ii])             #unique geneid vector
     unid<-geneid[ii]
     uniqueslide<-sld[ii]
     if(!is.null(weight))  Weight<-weight[ii]

   #reorder spike data
     id.s<-as.integer(unclass(factor(s.geneid)))
     slide.s<-as.integer(unclass(factor(s.sld)))
     ii<-order(slide.s,id.s)
     id.s<-id.s[ii]
     slide.s<-slide.s[ii]
     ratio.s<-s.rt[ii]
     intensity.s<-s.intn[ii]
     uniqueid.s<-unique(s.geneid[ii])
     unid.s<-s.geneid[ii]
     uniqueslide.s<-s.sld[ii]
     
     if(!is.null(s.weight)) Weight.s=s.weight[ii]
     
     geneno<-length(uniqueid)        #total gene number
     obsno<-length(id)               #number of row
     slidno<-length(unique(slide))   #slide number
     Nj<-as.vector(table(id))        #number of each clones in the experiment
     Ni<-as.vector(table(slide))     #number of clones in each slide
   #initial beta as the mean of ratio
    beta0<-beta<-.C("ymean",
        y=as.double(ratio),
      gno=as.integer(geneno),
       id=as.integer(id),
      obs=as.integer(obsno),
       nj=as.integer(Nj),
    ymean=double(geneno),PACKAGE="TwslmSpikeWeight")$ymean
    
    beta=beta0;beta.s0=0.0;
    Bfit<-NULL;Bfit.s<-NULL;
    fittedvalue<-rep(0,obsno)
    residi<-resid0i<-sqrt(sum(ratio^2)+sum(ratio.s^2))     #initialize resid0 and resid
    eps<-iter<-1

### sum weights for each gene for WLS regression
    if(!is.null(weight)){
      sumw<-.Call("sumweight",
                 as.double(Weight),
                 as.integer(id),
                 as.integer(geneno),PACKAGE="TwslmSpikeWeight")  #sum of weight for each gene
      sumw.s<-.Call("sumweight",
                 as.double(Weight.s),
                 as.integer(id.s),
                 as.integer(length(uniqueid.s)),PACKAGE="TwslmSpikeWeight")  #sum of weight for each spike
    }
  betas=median(ratio.s)  
  
while(eps>EPS){
  Bfit<-NULL;Bfit.s<-NULL;
  
  #B-spline by slide
  for(i in 1:length(unique(slide))){
     xstar=NULL;ystar=NULL
     betaresid<-.Call("resid", as.double(ratio[slide==i]-betas),
                               as.integer(id[slide==i]),
                               as.double(beta),PACKAGE="TwslmSpikeWeight")
     b<-intensity[slide==i]         ##non-spike 
     b.s<-intensity.s[slide.s==i]   ##spike
     bj<-bs(b,df=df,degree=degree)[1:length(b),1:df]
     bj.s<-bs(b.s,df=df,degree=degree)[1:length(b.s),1:df]
     bj<-cbind(rep(1,length(b)),bj)
     bj.s<-cbind(rep(1,length(b.s)),bj.s)
     
     if(is.null(weight)){
        xstar<-t(bj.s) %*% bj.s+t(bj)%*%bj
        ystar<-t(bj) %*% betaresid+t(bj.s) %*% (ratio.s[slide.s==i]-s.constant)
     }
     else {
       t1=t(bj * Weight[slide==i])
       t11=t1 %*% bj
       t2=t(bj.s*Weight.s[slide.s==i])
#xstar<-t(bj.s) %*% diag(Weight.s[slide.s==i]) %*% bj.s+t(bj) %*% diag(Weight[slide==i]) %*% bj
       xstar<-t2 %*% bj.s+t11
#ystar<-t(bj) %*% diag(Weight[slide==i]) %*% betaresid+t(bj.s) %*% diag(Weight.s[slide==i]) %*% (ratio.s[slide.s==i]-s.constant-betas)
       ystar<-t1 %*% betaresid+t2 %*% (ratio.s[slide.s==i]-s.constant)
      }
     xinverse<-solve(crossprod(xstar))
     lambdahat<-xinverse %*% t(xstar) %*% ystar
     Bfit<-c(Bfit,bj %*% lambdahat)
     Bfit.s<-c(Bfit.s, bj.s %*% lambdahat)

   }

   fittedvalue<-.Call("fittedvalue",as.double(Bfit),as.double(beta),as.integer(id),PACKAGE="TwslmSpikeWeight")

   residi<-sqrt(sum((ratio-fittedvalue)^2)+sum((ratio.s-Bfit.s-s.constant)^2))
   eps<-abs(residi-resid0i)/resid0i
cat(residi," ",resid0i," ",eps,"\n")   
   resid0i<-residi
   iter<-iter+1
   beta0<-beta


  #calculate beta value
 if(is.null(weight)){
    beta<-.C("ymean",
        y=as.double(ratio-Bfit-betas),
      gno=as.integer(geneno),
       id=as.integer(id),
      obs=as.integer(obsno),
       nj=as.integer(Nj),
      ymean=double(geneno),PACKAGE="TwslmSpikeWeight")$ymean
    betas<-median(ratio.s-Bfit.s,na.rm=TRUE)  
  }
  else {
   beta<-.Call("Wymean",
               as.double(Weight*(ratio-Bfit-betas)),
               as.integer(id),
               as.integer(geneno), 
               as.double(sumw),PACKAGE="TwslmSpikeWeight")
    
   beta.s<-.Call("Wymean",
                 as.double(Weight.s*(ratio.s-Bfit.s)),
                  as.integer(id.s),
                  as.integer(length(uniqueid.s)), 
                  as.double(sumw.s),PACKAGE="TwslmSpikeWeight")      
   betas<-median(beta.s,na.rm=TRUE)
  }
}
###***************  end of parameter estimation  **************

  beta0<-as.matrix(beta0)
  rownames(beta0)<-uniqueid
  colnames(beta0)<-"beta"

  parm<-list(name=uniqueid,beta=beta0,fittedvalue=fittedvalue, bfit=Bfit,
             slide=uniqueslide, id=unid,ratio=ratio, intensity=intensity,weight=Weight,
             spike.slide=uniqueslide.s,spike.id=unid.s,spike.bfit=Bfit.s,spike.ratio=ratio.s,
             spike.intensity=intensity.s,spike.weight=Weight.s,scale=NULL)

  cat("Calculation is end. ",date(),"\n")

return(parm)
}

#******************************************************************************
#                          Robust Model                                       *
#******************************************************************************

robust.spike <- function(sld,geneid,rt,intn,weight=NULL,s.sld,s.geneid,s.rt,s.intn,s.weight=NULL,s.constant=0.0,
      df=12,degree=3, robust.name="Tukey",scale.constant=2.5,weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)
      {

  cat("Calculation is running...",date(),"\n")
     EPS<-tol
     id<-as.integer(unclass(factor(geneid)))     #make id from 1,2,....J
     slide<-as.integer(unclass(factor(sld)))     #make slide 1,2,...,n
     ii<-order(slide,id)
     
     id<-id[ii]
     slide<-slide[ii]
     ratio<-rt[ii]
     intensity<-intn[ii]
     uniqueid<-unique(geneid[ii])             #unique geneid vector
     unid<-geneid[ii]
     uniqueslide<-sld[ii]
     if(!is.null(weight))  Weight<-weight[ii]
     else Weight<-rep(1,length(id))
     
   #reorder spike data
     id.s<-as.integer(unclass(factor(s.geneid)))
     slide.s<-as.integer(unclass(factor(s.sld)))
     ii<-order(slide.s,id.s)
     id.s<-id.s[ii]
     slide.s<-slide.s[ii]
     ratio.s<-s.rt[ii]
     intensity.s<-s.intn[ii]
     uniqueid.s<-unique(s.geneid[ii])
     unid.s<-s.geneid[ii]
     uniqueslide.s<-s.sld[ii]
     if(!is.null(s.weight)) Weight.s=s.weight[ii]
     else Weight.s=rep(1,length(id.s))
    
     geneno<-length(uniqueid)        #total gene number
     obsno<-length(id)               #number of row
     slidno<-length(unique(slide))   #slide number
     Nj<-as.vector(table(id))        #number of each clones in the experiment
     Ni<-as.vector(table(slide))     #number of clones in each slid
     beta<-beta0<-NULL               #initial beta as the mean of ratio
     if(is.null(ibeta)==TRUE){
       beta0<-.C("ymean",
              as.double(ratio),
              as.integer(geneno),
              as.integer(id),
              as.integer(obsno),
              as.integer(Nj),
              double(geneno),PACKAGE="TwslmSpikeWeight")[[6]]
      }
     else beta0<-ibeta  
   
    beta<-beta0
    Bfit<-NULL;Bfit.s<-NULL;
    fittedvalue<-rep(0,obsno)
    resid0i<-sqrt(sum(ratio^2)+sum(ratio.s^2))     #initialize resid0 and resid
    residi.s<-NULL;residi<-NULL;
    w<-Weight; ws<-Weight.s;   #initial weights as 1s
    eps<-1;
    sumw<-NULL
    isigma<-NULL 
    if(is.null(iscale)==TRUE) 
         isigma<-1.482602*median(abs(ratio-median(ratio)))               #initialize sigma
    else isigma<-iscale
    fsigma<-isigma
    betas=median(ratio.s)  

#****************while loop*****************************
   while(eps>EPS){

    Bfit<-NULL; Bfit.s<-NULL;
    for(i in unique(slide)){
        betaresid<-.Call("resid", as.double(ratio[slide==i]-betas),
                                  as.integer(id[slide==i]),
                                  as.double(beta),PACKAGE="TwslmSpikeWeight")

             b<-intensity[slide==i]
             b.s<-intensity.s[slide.s==i]
             bj<-bs(b,df=df,degree=degree)[1:length(b),1:df]
             bj.s<-bs(b.s,df=df,degree=degree)[1:length(b.s),1:df]         
             bj<-cbind(rep(1,length(b)),bj)
             bj.s<-cbind(rep(1,length(b.s)),bj.s)

             t1=t(bj * w[slide==i])
             t11=t1 %*% bj
             t2=t(bj.s*ws[slide.s==i])
             xstar<-t2 %*% bj.s+t11
             ystar<-t1 %*% betaresid+t2 %*% (ratio.s[slide.s==i]-s.constant)
             xinverse<-solve(crossprod(xstar))
             lambdahat<-xinverse %*% t(xstar) %*% ystar
             Bfit<-c(Bfit,bj %*% lambdahat)
             Bfit.s<-c(Bfit.s, bj.s %*% lambdahat)
    }

      fittedvalue <-.Call("fittedvalue",as.double(Bfit),as.double(beta),as.integer(id),PACKAGE="TwslmSpikeWeight")
      beta0<-beta
      residual<-ratio-fittedvalue             ##gene residual
      residual.s<-ratio.s-Bfit.s-s.constant   ##spike residual

      if(robust.name=="Huber"){
        residi<-.Call("Robustobj",
                      as.double(c(residual,residual.s)),
                      as.double(isigma),
                      as.double(weight.constant),
                      as.integer(1),PACKAGE="TwslmSpikeWeight")
        }                                #Huber's objective function
      else if(robust.name=="Tukey"){
        residi <- .Call("Robustobj",
                  as.double(c(residual,residual.s)),
                  as.double(isigma),
                  as.double(weight.constant),
                  as.integer(2),PACKAGE="TwslmSpikeWeight")
      }
      eps <- abs(residi-resid0i)/resid0i
      cat("eps=",eps,"resid0i=",resid0i,"residi=",residi,"\n")
      resid0i<-residi
      isigma<-fsigma               #assign fisigma to isigma
     
    #Robust regression
        if(robust.name=="Tukey"){
          fsigma<-.Call("Tukeyscale",
                    as.double(c(residual,residual.s)),
                        as.integer((df+1)*slidno+geneno),
                        as.double(isigma),
                        as.double(scale.constant),PACKAGE="TwslmSpikeWeight")
          w<-.Call("TukeyWeight",as.double(residual),as.double(fsigma),as.double(weight.constant),PACKAGE="TwslmSpikeWeight")
          w<-w*Weight
          ws<-.Call("TukeyWeight",as.double(residual.s),as.double(fsigma),as.double(weight.constant),PACKAGE="TwslmSpikeWeight")
          ws<-ws*Weight.s
         }
         else if (robust.name=="Huber"){
            fsigma<-.Call("Huberscale",
                     as.double(c(residual,residual.s)),
                         as.integer((df+1)*slidno+geneno),
                         as.double(isigma),
                         as.double(scale.constant),PACKAGE="TwslmSpikeWeight")
             fsigma<-sqrt(fsigma)
             w<-.C("Huber",
                 as.double(residual),
                 double(obsno),
                 as.integer(obsno),
                 as.double(fsigma),
                 as.double(weight.constant),PACKAGE="TwslmSpikeWeight")[[2]]
             w<-w*Weight
             ws<-.C("Huber",
                 as.double(residual.s),
                 double(length(residual.s)),
                 as.integer(length(residual.s)),
                 as.double(fsigma),
                 as.double(weight.constant),PACKAGE="TwslmSpikeWeight")[[2]] 
             ws<-ws*Weight.s    
        }

      sumw<-.Call("sumweight",
                 as.double(w),
                 as.integer(id),
                 as.integer(geneno),PACKAGE="TwslmSpikeWeight")  #sum of weight for each gene
      sumws<-.Call("sumweight",
                 as.double(ws),
                 as.integer(id.s),
                 as.integer(length(uniqueid.s)),PACKAGE="TwslmSpikeWeight")  #sum of weight for each gene
      beta<-.Call("Wymean",
                 as.double(w*(ratio-Bfit-betas)),
                 as.integer(id),
                 as.integer(geneno), 
                 as.double(sumw),PACKAGE="TwslmSpikeWeight")
      beta.s<-.Call("Wymean",
                as.double(ws*(ratio.s-Bfit.s)),
                 as.integer(id.s),
                 as.integer(length(uniqueid.s)), 
                 as.double(sumws),PACKAGE="TwslmSpikeWeight")      
      betas<-median(beta.s,na.rm=TRUE)
}
  
#save the results in
  beta0<-as.matrix(beta0)
  rownames(beta0)<-uniqueid
  colnames(beta0)<-"beta"

    parm<-list(name=uniqueid,beta=beta0,fittedvalue=fittedvalue, bfit=Bfit,
          slide=uniqueslide, id=unid,ratio=ratio,intensity=intensity, weight=Weight,rweight=w/Weight,
          spike.slide=uniqueslide.s,spike.id=unid.s,spike.bfit=Bfit.s,spike.ratio=ratio.s,
          spike.intensity=intensity.s,spike.weight=Weight.s,rspike.weight=ws/Weight.s,scale=isigma^2)

  cat("Calculation is end. ",date(),"\n")    
  return(parm)
}

#******************************************************
#   Blockwise Normalization Function                  *
#******************************************************

BlockByBlock.spike <- function(sld,blk,geneid,rt,intn,weight=NULL,
   s.sld,s.blk,s.geneid,s.rt,s.intn,s.weight=NULL,s.constant=0.0,
   df=12,degree=3,robust=TRUE,robust.name="Tukey",scale.constant=2.5, 
   weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5){
  
 beta<-NULL;slide<-NULL; block<-NULL;ID<-NULL;
 ratio<-NULL;intensity<-NULL; bfit<-NULL; fittedvalue<-NULL;scale<-NULL;name<-NULL;
 spike.slide<-NULL;spike.block<-NULL;spike.ratio<-NULL;spike.intensity<-NULL;spike.id<-NULL;
 spike.bfit<-NULL;Weight=NULL;spike.weight=NULL;rweight=NULL;rspike.weight=NULL
 
 blockid<-unique(cbind(blk,geneid))
 cat("Block Number=")
 for(i in unique(blk)){
   cat(i," ")
   as<-sld[blk==i]
   ai<-geneid[blk==i]
   ar<-rt[blk==i]
   ain<-intn[blk==i]
   aw<-weight[blk==i]
   bi<-rep(i,length(as))
 ##spikes
   sas<-s.sld[s.blk==i]
   sai<-s.geneid[s.blk==i]
   sar<-s.rt[s.blk==i]
   sain<-s.intn[s.blk==i]
   saw<-s.weight[s.blk==i]
   sbi<-rep(i,length(sas))
   
   a<-NULL
   ii<-(blockid[,1]==i)
   if(robust==TRUE){
     a<-robust.spike(sld=as,geneid=ai,rt=ar,intn=ain,weight=aw,
          s.sld=sas,s.geneid=sai,s.rt=sar,s.intn=sain,s.weight=saw,s.constant=s.constant,
          df=df,degree=degree,robust.name=robust.name, scale.constant=scale.constant,
          weight.constant=weight.constant,ibeta=ibeta[ii],iscale=iscale[i], tol=tol)
   }   
   else
   {
    a<-spike(sld=as,geneid=ai,rt=ar,intn=ain,weight=aw,
      s.sld=sas,s.geneid=sai,s.rt=sar,s.intn=sain,s.weight=saw,
      s.constant=s.constant,df=df,degree=degree,tol=tol)
   }  
   beta<-c(beta,a$beta)
   
   slide<-c(slide,a$slide)
   block<-c(block,bi)
   ID<-c(ID,a$id)
   ratio<-c(ratio,a$ratio)
   intensity<-c(intensity,a$intensity)
   bfit<-c(bfit,a$bfit)
   fittedvalue<-c(fittedvalue,a$fittedvalue)
   Weight=c(Weight,a$weight)
   rweight=c(rweight,a$rweight)
   
   spike.block<-c(spike.block,sbi)
   spike.id<-c(spike.id,a$spike.id)
   spike.ratio<-c(spike.ratio,a$spike.ratio)
   spike.intensity<-c(spike.intensity,a$spike.intensity)
   spike.bfit<-c(spike.bfit,a$spike.bfit) 
   spike.slide<-c(spike.slide,a$spike.slide)
   spike.weight<-c(spike.weight,a$spike.weight)
   rspike.weight<-c(rspike.weight,a$rspike.weight)
   
   scale<-c(scale,a$scale)
   name<-c(name,a$name)
 }
   beta<-as.matrix(beta)
   rownames(beta)<-name
   colnames(beta)<-"beta"
 
 result<-list(name=name,beta=beta,bfit=bfit,fittedvalues=fittedvalue,slide=slide,id=ID,block=block,
             ratio=ratio,intensity=intensity,weight=Weight,rweight=rweight,spike.slide=spike.slide,spike.id=spike.id,
             spike.block=spike.block,spike.ratio=spike.ratio,spike.intensity=spike.intensity,
             spike.bfit=spike.bfit,spike.weight=spike.weight,rspike.weight=rspike.weight,scale=scale)

 return(result)
}


#***************************************************
# Semiparametric model normalization main function *
#***************************************************

TwslmSpikeWeight <- function(sld,blk,geneid,rt,intn,weight,s.sld,s.blk,s.geneid,s.rt,s.intn,s.weight,
             s.constant=0.0,df=12,degree=3,block.norm=FALSE,robust=TRUE,robust.name="Tukey",
             scale.constant=2.5,weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)
{
# cat("Start", date()," ")
  param<-NULL
  if(block.norm==FALSE){
    if(robust==FALSE)
        param <- spike(sld=sld,geneid=geneid,rt=rt,intn=intn,weight=weight,
        s.sld=s.sld,s.geneid=s.geneid,s.rt=s.rt,s.intn=s.intn,s.weight=s.weight,
        s.constant=s.constant,df=df,degree=degree,tol=tol)
    else {
    param <- robust.spike(sld=sld,geneid=geneid,rt=rt,intn=intn,weight=weight,
    s.sld=s.sld,s.geneid=s.geneid,s.rt=s.rt,s.intn=s.intn,s.weight=s.weight,
    s.constant=s.constant,df=df,degree=degree,robust.name=robust.name,
    scale.constant=scale.constant,weight.constant=weight.constant,ibeta=ibeta,iscale=iscale,tol=tol)
   }
  }
  else param <- BlockByBlock.spike(sld=sld,blk=blk,geneid=geneid,rt=rt,intn=intn,weight=weight,
     s.sld=s.sld,s.blk=s.blk,s.geneid=s.geneid,s.rt=s.rt,s.intn=s.intn,s.weight=s.weight,
     s.constant=s.constant,df=df,degree=degree,robust=robust,robust.name=robust.name,
     scale.constant=scale.constant,weight.constant=weight.constant,ibeta=ibeta,iscale=iscale,tol=tol)
 
 cat("End",date(),"\n")
 
 return(param)
}
