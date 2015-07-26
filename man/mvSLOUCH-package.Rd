\name{mvSLOUCH-package}
\alias{mvSLOUCH-package}
\alias{mvSLOUCH}
\docType{package}
\title{Multivariate Ornstein-Uhlenbeck type stochastic differential equation models 
for phylogenetic comparative data.
}
\description{
The package allows for  maximum likelihood estimation, simulation and study of properties 
of multivariate Brownian motion \deqn{\begin{array}{rcl}dX(t) & = & \Sigma dB(t),\end{array}}{dX(t) = S dB(t),}
OU \deqn{\begin{array}{rcl}dY(t) & = & -A(Y(t)-\Psi(t))dt + \Sigma dB(t)\end{array}}{dY(t) = -A(Y(t)-Psi(t))dt + SdB(t)}
and OUBM  
\deqn{\begin{array}{rcl}
dY(t) & = & -A(Y(t)-\Psi(t)- A^{-1}BX(t))dt + \Sigma_{yy} dB(t) \\ dX(t) & = & \Sigma_{xx} dB(t)
\end{array}
}{
dY(t)  =  -A(Y(t)-Psi(t)- solve(A)BX(t))dt + Syy dB(t) 

dX(t)  = Sxx dB(t)
}
models that evolve on a phylogenetic tree.

This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
Please understand that there may still be bugs and errors. Use it at your own risk. 
We take no responsibility for any errors or omissions in this package or for any misfortune 
that may befall you or others as a result of its use. Please send comments and report 
bugs to Krzysztof Bartoszek at bartoszekkj@gmail.com .
}
\details{
\tabular{ll}{
Package: \tab mvSLOUCH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-26\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The package allows for  maximum likelihood estimation, simulation and study of properties 
of multivariate Brownian motion 
\deqn{\begin{array}{rcl}dX(t) & = & \Sigma dB(t),\end{array}}{dX(t) = S dB(t),}
OU \deqn{\begin{array}{rcl}dY(t) & = & -A(Y(t)-\Psi(t))dt + \Sigma dB(t)\end{array}}{dY(t) = -A(Y(t)-Psi(t))dt + SdB(t)}
and OUBM  
\deqn{\begin{array}{rcl}
dY(t) & = & -A(Y(t)-\Psi(t)- A^{-1}BX(t))dt + \Sigma_{yy} dB(t) \\ dX(t) & = & \Sigma_{xx} dB(t)
\end{array}
}{
dY(t)  =  -A(Y(t)-Psi(t)- solve(A)BX(t))dt + Syy dB(t) 

dX(t)  = Sxx dB(t)
}
models that evolve on a phylogenetic tree.

The estimation functions are \code{BrownianMotionModel}, \code{ouchModel} (OUOU) and 
\code{mvslouchModel} (mvOUBM).
The rely on a combination of least squares and numerical optimization techniques. 
A wrapper funciotn for all of them is \code{estimate.evolutionary.model}, in tries all
three model with different matrix parameter classes and then returns the best model 
based on the AICc.

The simulation functions are simulBMProcPhylTree, simulOUCHProcPhylTree, 
simulMVSLOUCHProcPhylTree.

The phylogeny provided to them is required to be of the same type that the \code{ouch} 
package requires.

The package uses the functions .sym.par and .sym.unpar from the \code{ouch} package to 
parametrize symmetric matrices. The package also uses the internal structure of 
\code{ouch}'s (as of version 2.8.2) tree object.

In the case the mvOUBM model with a single response trait the package \code{slouch} is a 
recommended alternative.


}

\author{
Krzysztof Bartoszek
Maintainer: <bartoszekkj@gmail.com>
}
\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.

Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

Felsenstein, J. (1985) 
  Phylogenies and the comparative method.
  American Naturalist 125:1-15.

Hansen, T.F. (1997) 
  Stabilizing selection and the comparative analysis of adaptation. 
  Evolution 51:1341-1351.

Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and biological 
  errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.

Hansen, T.F. and Pienaar, J. and Orzack, S.H. (2008) 
  A comparative method for studying adaptation to randomly evolving environment. 
  Evolution 62:1965-1977.

Labra, A., Pienaar, J. & Hansen, T.F. (2009) 
  Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and 
  niche tracking. 
  The American Naturalist 174:204-220.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}

\keyword{models}

\examples{
\dontrun{ ##It takes too long to run this
### We will first simulate a small phylogenetic tree using functions from ape and ouch.
### For simulating the tree one could also use alternative functions, 
## eg. sim.bd.taxa from the TreeSim package
phyltree<-ape2ouch(rtree(5))

### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
1:(phyltree@nnodes-phyltree@nterm))

### Define a vector of regimes.
regimes<-c("small","small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the different models.
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),
"large"=c(-1,1,0.5)),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))
OUBMparameters<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,0),c(0,5)),
B=matrix(c(2,-2),ncol=1,nrow=2),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(0,1,1),Sxx=matrix(1,1,1),
Syx=matrix(0,ncol=1,nrow=2),Sxy=matrix(0,ncol=2,nrow=1))

### Now simulate the data and remove the values corresponding to the internal nodes.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[-(1:(phyltree@nnodes-phyltree@nterm)),]
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[-(1:(phyltree@nnodes-phyltree@nterm)),]
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### Recover the parameters of the SDEs.
BMestim<-BrownianMotionModel(phyltree,BMdata)
OUOUestim<-ouchModel(phyltree,OUOUdata,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")

### Use the wrapper function
estimResultsBM<-estimate.evolutionary.model(demotreeouch,BMdata,regimes=NULL,
root.regime=NULL,M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),kY=2,
doPrint=TRUE)
estimResultsOUOU<-estimate.evolutionary.model(demotreeouch,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),kY=2,
doPrint=TRUE)
estimResultsOUBM<-estimate.evolutionary.model(demotreeouch,OUBMdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),kY=2,
doPrint=TRUE)
## In the wrapper function the resulting best found model parameters are in
## estimResultsBM$BestModel$ParamsInModel
## estimResultsOUOU$BestModel$ParamsInModel
## estimResultsOUBM$BestModel$ParamsInModel

### And summarize them.
BM.summary<-SummarizeBM(phyltree,BMdata,BMestim$ParamsInModel,t=c(1),
dof=BMestim$ParamSummary$dof,calcCI=FALSE)
OUOU.summary<-SummarizeOUCH(phyltree,OUOUdata,OUOUestim$FinalFound$ParamsInModel,
regimes,t=c(1),dof=OUOUestim$FinalFound$ParamSummary$dof,calcCI=FALSE)
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMestim$FinalFound$ParamsInModel,
regimes,t=c(1),dof=OUBMestim$FinalFound$ParamSummary$dof,calcCI=FALSE)
### if one would want the confidence intervals then set calcCI=TRUE
}
}
