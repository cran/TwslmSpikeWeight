\name{TwslmSpikeWeight}
\alias{TwslmSpikeWeight}
\alias{spike}
\alias{robust.spike}
\alias{BlockByBlock.spike}


\title{Normalization of cDNA microarray data using spike and spot quality 
information in the two-way semi-linear model (TW-SLM) }
\description{
By incoporating spike and spot quality information, normalization of
cDNA microarray data using the TW-SLM could be improved. Two methods
are available for estimation, including robust estimation methods and the 
least square method. The B-splines is used to estimate nonparametric normalization 
curves in the model.
}
\usage{
TwslmSpikeWeight(sld,blk,geneid,rt,intn,weight,s.sld,s.blk,s.geneid,s.rt,
            s.intn,s.weight,s.constant=0.0,df=12,degree=3,block.norm=FALSE,
            robust=TRUE,robust.name="Tukey", scale.constant=2.5,
            weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)

spike(sld,geneid,rt,intn,weight=NULL,s.sld,s.geneid,s.rt,s.intn,
      s.weight=NULL,s.constant=0.0,df=12,degree=3,tol=1e-5)
                
robust.spike(sld,geneid,rt,intn,weight=NULL,s.sld,s.geneid,s.rt,s.intn,
             s.weight=NULL,s.constant=0.0,df=12,degree=3,
             robust.name="Tukey",scale.constant=2.5,weight.constant=4.685,
             ibeta=NULL,iscale=NULL,tol=1e-5)
                
BlockByBlock.spike(sld,blk,geneid,rt,intn,weight=NULL,s.sld,s.blk,s.geneid,
                   s.rt,s.intn,s.weight=NULL,s.constant=0.0,df=12,degree=3,
                   robust=TRUE,robust.name="Tukey",scale.constant=2.5, 
                   weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=1e-5)
}
\arguments{
\item{sld}{
a vector of array or slide numbers. This argument is required.
}
\item{blk}{
a vector of block numbers. This argument is required only for blockwise normalization.
}
\item{geneid}{
a vector of gene identification numbers, can be numerical numbers or gene names.
This argument is required.
}
\item{rt}{
a vector of \eqn{log_2} intensity ratio, i.e. \eqn{log_2(Cy5/Cy3)}. This argument is required.
}
\item{intn}{
a vector of average of log two total intensity, i.e. \eqn{0.5log_2(Cy5*Cy3)}. 
This argument is required.
}
\item{weight}{a vector of weights assigned for non-spike spots. This argument 
is required for normalization using spot quality information or pre-assigned 
weights for non-spikes spots are provided. Weights must be greater or equal to zero.
}
\item{s.sld}{a vector of array or slide numbers for spike spots. This argument is required.
}
\item{s.blk}{
a vector of block numbers where spike spots locate. This argument is required only for 
blockwise normalization.
}
\item{s.geneid}{
a vector of spike identification numbers, can be numerical numbers or spike names.
This argument is required.
}
\item{s.rt}{
a vector of \eqn{log_2} intensity ratio for spikes, i.e. \eqn{log_2(Cy5/Cy3)}. 
This argument is required.
}
\item{s.intn}{
a vector of averages of log two total intensity for spikes, i.e. \eqn{0.5log_2(Cy5*Cy3)}. 
This argument is required.
}
\item{s.weight}{
a vector of weights assigned for spike spots. This argument is required for 
normalization using spot quality information or pre-assigned 
weights for spikes spots are provided. Weights must be greater or equal to zero.
}
\item{s.constant}{
a constant which spikes centered after normalization. The default value is zero. 
This argument is required.
}
\item{df}{
the degrees of freedom for B-spline smooth. The default is 12.
}
\item{degree}{
the order of polynomials in the B-splines. The default is 3, the cubic spline.
}
\item{block.norm}{
a logical value indicating whether blockwise normalization is performed or not.
The default is FALSE, which means the default normalization is slide by slide normalization.
}
\item{robust}{
a logical value indicating if the robust procedure is incoporated in normalization.
The default is TRUE, which means normalization is conducted using a robust method
in the TW-SLM. The least square method is used if this argument is FALSE.
}
\item{robust.name}{
a name for the robust procedure. The default is "Tukey", which means the
location and scale parameters are estimated iteratively with Tukey's bisquare
weight function. Another option is "Huber", which uses Huber's
weight function and the location and scale parameters are estimated iteratively.
This option works only if the robust argument is TRUE.
}
\item{scale.constant}{
a constant chosen for scale estimation in the robust TW-SLM. The default is 2.5.
}
\item{weight.constant}{
a constant chosen for robust location estimation. The default is 1.345 for Huber's 
weight function and 4.685 for Tukey's weight function.
}
\item{ibeta}{
a vector for initalization of \eqn{\beta}. The default is NULL. The ordinary least square 
estimators for \eqn{\beta} is one choice of starting values for the robust TW-SLM. 
Giving this value will speed up convergence.
}
\item{iscale}{
a value for initalization of the scale parameter in the robust model. The default is NULL.
Giving this value will speed up convergence of the algorithm.
}
\item{tol}{
a convergent criteria for iterative estimation procedure. The default is 1e-5.
}
}
\details{
Normalization is a basic step in the analysis of microarray data.
Widely used normalization method for cDNA microarray data was the {\it loess}
normalization method proposed by Yang et al.(2001). This method requires 
that at least one of the two underlying biological assumptions, i.e. either
(i) a small fraction of genes in the experiment are differentially expressed; or
(ii) the up-regulated genes and the down-regulated genes are distributed symmetrically.
The TW-SLM is a generalization of the semiparametric 
regression model. It does not require either of the above two assumptions 
for normalization of cDNA microarray data. Side information including spike information
and spot quality information could be used to improve normalization of cDNA microarray
data using TW-SLM.

The TW-SLM has the form
\deqn{y_{ij(s)}=\phi_i(x_{ij(s)})+\beta_j(s)+\epsilon_{ij(s)}.}
where \eqn{y_{ij(s)}=log_2(Cy5/Cy3)}, \eqn{\phi_i(x_{ij(s)})} is the intensity
dependent normalization curve for slide \eqn{i}, \eqn{x_{ij(s)}=0.5log_2(Cy5*Cy3)},
\eqn{\beta_j} is the relative effect of gene \eqn{j}, \eqn{\epsilon_{ij(s)}} is
the residual term, for \eqn{i=1,\ldots,n}, where \eqn{n} is the total number of slides,
\eqn{j=1,\ldots,J}, where \eqn{J} is the total number of genes in the experiment.
Applying the TW-SLM to spike spots specifically and integrate spikes and non-spike genes
together in the normalization process. Spikes and non-spike genes share common 
normalization curves. Spot quality information could also be incorporated in 
the TW-SLM so that normalization might be improved. Details of integration of 
spikes and spot quality information in normalization of cDNA microarray data using
the TW-SLM can be found in Wang et al. (2007).

The \code{TwslmSpikeWeight} package implements the TW-SLM for normalization
of cDNA microarray data considering spike spots information and/or spot quality
 information. Two robust estimation procedures are implemented in the 
 current version of \code{TwslmSpikeWeight}: Huber's method (1981) and Tukey's 
 method (1986). The \code{spike} function implements the ordinary least square
 estimation method and the \code{robust.spike} function implements robust 
 estimation methods. The \code{BlockByBlock.spike} function carries out
 block-wise normalization.
}
\value{
  An object of a list is returned with components:
\item{name}{
a vector of names of unique genes.
}
\item{beta}{
an estimated parameters of relative gene expression level for each gene.
}
\item{fittedvalue}{a vector of fitted values in the TW-SLM.
}
\item{bfit}{
a vector of fitted values for normalization curves.
}
\item{slide}{
a vector of slide number from the input the function. The order is
different from the input "sld" vector.
}
\item{id}{
a vector of gene ID from the input vector "geneid" with a different order.
}
\item{ratio}{
a vector of the log two intensity ratio from the input vector "rt" with a different
order.
}
\item{intensity}{
a vector of average log two total intensity from the input vector "intn" with
a different order.
}
\item{weight}{
a vector of weights used in the TW-SLM. These weights could be 
pre-assigned weights or calculated weights using spot quality information.
}
\item{rweight}{
a vector of weights calculated at the last step of convergence in  
the robust TW-SLM.
}
\item{spike.slide}{
a vector of slide number for spikes. The order is
different from the input vector "s.sld".
}
\item{spike.id}{
a vector of spike identification numbers from the input vector "s.geneid" 
with a different order.
}
\item{spike.ratio}{
a vector of the log two intensity ratio for spikes' intensity with 
a different order from input "s.rt".
}
\item{spike.intensity}{
a vector of average log two total intensity for spikes with a different 
order from input vector "s.intn".
}
\item{spike.bfit}{
a vector of fitted values for spikes after normalization.
}
\item{spike.weight}{
a vector of weights for spikes in the TW-SLM. These weights could be 
pre-assigned weights or calculated weights using spike spot quality information.
}
\item{rspike.weight}{
a vector of weights for spikes calculated at the last step of convergence in the 
robust TW-SLM.
}
\item{scale}{
a scale estimator in the two-way semilinear model.
}
}
\note{ \code{TwslmSpikeWeight} is the main function to control which normalizatin method will 
be used. \code{spike} is the function for the TW-SLM using the ordinary
least squares, \code{robust.spike} is the function for robust estimation of the TW-SLM,
\code{BlockByBlock.spike} is the function for blockwise normalization using spikes.
}

\references{
   Huang, J., Wang, D. & Zhang, C.H. (2005), 
   \bold{A Two-way Semi-Linear Model for Normalization and Analysis of Microarray Data}.
\emph{Journal of the American Statistical Association, 100(471):814-829}

   Wang, D., Huang, J., Xie, H., Manzella, L., Soares, M. B.,  \bold{A robust two-way semi-linear model for
 normalization of cDNA microarray data},\emph{BMC Bioinformatics 2005, 6:14}.

  Wang, D., Zhang, C-H,Soares, M. B., Huang, J.(2007) \bold{Systematic approaches for 
  incorporating control spots and data quality information to improve 
  normalization of cDNA microarray data}, \emph{Journal of Biopharmaceutical 
  Statistics, 17:415-431}.
  
 Yang, Y. H., Dudoit, S., Luu, P. & Speed, T. P. (2001), \bold{Normalization for cDNA microarray}.
   In Bittner, M.L., Chen, Y., Dorsel, A.N. and Dougherty, E.R.(eds), Microarrays: Optical Technologies and
   Informatics. SPIE, Society for Optical Engineering, San Jose, CA, 4266.

   Huber, P.J. (1981), \bold{Robust Statistics}, John Wiley & Sons.

   Hampel, F. R., Ronchetti, E. M., Rousseeuw, P. J. & Stahel, W. A.(1986), \bold{Robust Statistics-The Approach
   Based on Influence Functions}, John Wiley & Sons.
}
\examples{

## The example used in the paper Wang et al. (2007) will be used to 
## demonstrate utilization of the developed package for normlization of cDNA 
## microarray data by incorporating spike and spot quality information.

## load data used for the calculation

#require(TwslmSpikeWeight)

data(HU,"HUweightdata")

## a function to calculate weights of spots

ratioweight <- function(RF,RB,SDRF,SDRB,GF,GB,SDGF,SDGB,NF,NB){

  weight <- 1/((SDRF^2/NF+SDRB^2/NB)/(RF-RB)^2+(SDGF^2/NF+SDGB^2/NB)/(GF-GB)^2)

  return(weight)
}

## spikes ids
zerospike <- c("01-NBB-f-02","02-NAP-a-01","03-NAP-a-05","04-NAP-a-09","05-NXEg-e-10",
"06-NXEg-h-07","07-NBB-h-08","08-NBB-f-05","09-NAP-b-02","10-NS8W-d-02","11-NS8W-c-10",
"12-NAP-c-05","13-NAP-b-06","14-NAP-b-08","15-NAP-d-08","16-NAP-c-02","17-NBB-g-12",
"18-NAP-d-03","19-NAP-d-06","20-NXEg-h-06","21-NAP-c-08","22-NBB-h-07","23-NBB-e-09",
"24-NAP-d-09","25-NAP-b-09","26-NBB-f-08","SC1","SC3","SC5")


## calculate weights

allw <- ratioweight(RF=HUweightdata$RF,RB=HUweightdata$RB,SDRF=HUweightdata$SDRF,
  SDRB=HUweightdata$SDRB,GF=HUweightdata$GF,GB=HUweightdata$GB,SDGF=HUweightdata$SDGF,
  SDGB=HUweightdata$SDGB,NF=HUweightdata$NF,NB=HUweightdata$NB)

HUid <- HU$HUid
HUratio <- HU$HUratio
HUint <- HU$HUint
HUslide <- HU$HUslide
HUblock <- HU$HUblock

flag <- match(HUid,zerospike,nomatch=0)>0

## calculate weights

w <- allw[!flag]
huid <- HUid[!flag]
hurt <- HUratio[!flag]
huint <- HUint[!flag]
husld <- HUslide[!flag]
hublk <- HUblock[!flag]

##spikes information

sw <- allw[flag]
shuid <- HUid[flag]
shurt <- HUratio[flag]
shuint <- HUint[flag]
shusld <- HUslide[flag]
shublk <- HUblock[flag]

## spike,robust,weight;

huSRW <- TwslmSpikeWeight(sld=husld,geneid=huid,rt=hurt,intn=huint, weight=w,
 s.sld=shusld,s.geneid=shuid,s.rt=shurt,s.intn=shuint,s.weight=sw,s.constant=0.0,
 df=12,degree=3,block.norm=FALSE,robust=TRUE,robust.name="Tukey", scale.constant=2.5,
 weight.constant=4.685,ibeta=NULL,iscale=NULL,tol=0.0001)

a <- huSRW$slide
b <- huSRW$spike.slide
rr <- huSRW$ratio[a==1]
intt <- huSRW$intensity[a==1]
ff <- huSRW$bfit[a==1]
ii <- order(intt)

plot(intt,rr,xlim=c(0,16))
lines(intt[ii],ff[ii],col="red")
points(huSRW$spike.intensity[b==1],huSRW$spike.ratio[b==1],col="blue")

##zerospike estimates( mean median, sd)

zerospikemeanR <- aggregate(huSRW$spike.ratio-huSRW$spike.bfit,by=list(huSRW$spike.id),mean)
spikemeanR <- mean(2^zerospikemeanR$x)
spikemedianR <- median(2^zerospikemeanR$x)
spikesdR <- sd(2^zerospikemeanR$x)

}
\keyword{robust}
\keyword{smooth}
\keyword{nonparametric}
\author{
    Deli Wang \email{deli.wang@ccc.uab.edu} Jian Huang \email{jian@stat.uiowa.edu}
}


