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
dY(t) & = & -A(Y(t)-(\Psi(t)- A^{-1}BX(t)))dt + \Sigma_{yy} dB(t) \\ dX(t) & = & \Sigma_{xx} dB(t)
\end{array}
}{
dY(t)  =  -A(Y(t)-(Psi(t)- solve(A)BX(t)))dt + Syy dB(t) 

dX(t)  = Sxx dB(t)
}
models that evolve on a phylogenetic tree.

This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
Please understand that there may still be bugs and errors. Use it at your own risk. 
We take no responsibility for any errors or omissions in this package or for any misfortune 
that may befall you or others as a result of its use. Please send comments and report 
bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
}
\details{
\tabular{ll}{
Package: \tab mvSLOUCH\cr
Type: \tab Package\cr
Version: \tab 2.7.6\cr
Date: \tab 2023-11-19\cr
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
\code{mvslouchModel} (mvOUBM). They rely on a combination of least squares and numerical 
optimization techniques. A wrapper function for all of them is \code{estimate.evolutionary.model}, 
it tries all three models with different matrix parameter classes and then returns the best model 
based on the AICc.

The simulation functions are simulBMProcPhylTree, simulOUCHProcPhylTree, 
simulMVSLOUCHProcPhylTree.

The phylogeny provided to them should be of the \code{phylo} (package \pkg{ape}) format.

The package uses the functions \code{.sym.par()} and \code{.sym.unpar()} from the 
\pkg{ouch} package to parametrize symmetric matrices. 

In the case the mvOUBM model with a single response trait the package \pkg{slouch} is a 
recommended alternative.

The package uses \pkg{PCMBase}'s \code{PCMLik()} function as the engine to do calculate
the likelihood and phylogenetic least squares. If the \pkg{PCMBaseCpp} package is 
installed \pkg{mvSLOUCH} can take advantage of it to significantly decrease
the running time. The \pkg{PCMBaseCpp} package is available from 
\url{https://github.com/venelin/PCMBaseCpp}.
}

\author{
Krzysztof Bartoszek
Maintainer: <krzbar@protonmail.ch>
}
\references{
Bartoszek, K. and Fuentes-Gonzalez, J. and Mitov, V. and Pienaar, J. and Piwczynski, M.
  and Puchalka, R. and Spalik, K. and Voje, K. L. (2023) 
  Model Selection Performance in Phylogenetic Comparative Methods Under Multivariate 
  Ornstein-Uhlenbeck Models of Trait Evolution, Systematic Biology 72(2):275-293.

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

Mitov, V. and Bartoszek, K. and Asimomitis, G. and Stadler, T. (2020) 
  Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts
  Theoretical Population Biology 131:66-78.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}

\keyword{models}
\keyword{htest}
\keyword{multivariate}
\keyword{manip}
\keyword{hplot}
\keyword{datagen}
\keyword{array}
\keyword{algebra}

\examples{
## Set the version of the random number generator
## so that the results are always reproducible
RNGversion(min(as.character(getRversion()),"3.6.1"))
## Set the random seed, and random number generator
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, 
## e.g. sim.bd.taxa from the TreeSim package
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed.
## It enhances the phylo object with multiple additional fields, e.g.,
## the lineage for each tip, that are used by subsequent mvSLOUCH functions.
phyltree<-phyltree_paths(phyltree)

## Define a vector of regimes. There are two of them: small and large;
## their layout is random.
regimes<-c("small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the different models.
### There is no particular meaning attached to these parameters.
### BM parameters
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))
### OUOU parameters
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),
"large"=c(-1,1,0.5)),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))
### OUBM parameters
OUBMparameters<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,0),c(0,5)),
B=matrix(c(2,-2),ncol=1,nrow=2),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(0,1,1),Sxx=matrix(1,1,1),
Syx=matrix(0,ncol=1,nrow=2),Sxy=matrix(0,ncol=2,nrow=1))

### Now simulate the data under 3D BM, OUOU and OUBM (the third variable is the
### BM one) models of evolution.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
### extract just the values at the tips of the phylogeny
BMdata<-BMdata[phyltree$tip.label,,drop=FALSE]
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
### extract just the values at the tips of the phylogeny
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
### extract just the values at the tips of the phylogeny
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

### Recover the parameters of the SDEs.
### Recover under the BM model.
BMestim<-BrownianMotionModel(phyltree,BMdata)
### reset to the original version of the random number generator.
RNGversion(as.character(getRversion()))
\dontrun{ 
### It takes too long to run this from this point.
### We do not reduce the number of iterations of the optimier
### (as we did in other manual page examples, so the running
### times will be as with default settings.
### Recover under the OUOU model.
OUOUestim<-ouchModel(phyltree,OUOUdata,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")
### Recover under the OUBM model.
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")

### Usage of the wrapper function that allows for estimation under multiple 
### models in one call, here we use the default collection of models offered 
### by mvSLOUCH. This includes BM, a number of OUOU and OUBM models.
### For data simulated under BM.
estimResultsBM<-estimate.evolutionary.model(phyltree,BMdata,regimes=NULL,
root.regime=NULL,M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),
kY=2,doPrint=TRUE)
### For data simulated under OUOU.
estimResultsOUOU<-estimate.evolutionary.model(phyltree,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),
kY=2,doPrint=TRUE)
### For data simulated under OUBM.
estimResultsOUBM<-estimate.evolutionary.model(phyltree,OUBMdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),
kY=2,doPrint=TRUE)
## In the wrapper function the resulting best found model parameters are in
## estimResultsBM$BestModel$ParamsInModel
## estimResultsOUOU$BestModel$ParamsInModel
## estimResultsOUBM$BestModel$ParamsInModel

### Summarize the best found models.
### Under the BM model for BM data. One needs to check whether the 
### model in BMestim$ParamsInModel corresponds to a BM model (here it does).
BM.summary<-SummarizeBM(phyltree,BMdata,BMestim$ParamsInModel,t=c(1),
dof=BMestim$ParamSummary$dof)
### Under the OUOU model for OUOU data. One needs to check whether the 
### model in OUOUestim$ParamsInModel corresponds to an OUOU model (here it does).
OUOU.summary<-SummarizeOUCH(phyltree,OUOUdata,OUOUestim$FinalFound$ParamsInModel,
regimes,t=c(1),dof=OUOUestim$FinalFound$ParamSummary$dof)
### Under the OUBM model for OUBM data. One needs to check whether the 
### model in OUBMestim$ParamsInModel corresponds to a OUBM model (here it does).
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMestim$FinalFound$ParamsInModel,
regimes,t=c(1),dof=OUBMestim$FinalFound$ParamSummary$dof)


### Now run the parametric bootstrap to obtain confidence intervals for some parameters.
### For the BM model, we want CIs for the ancestral state, and the "sqaure" of the 
### diffusion (evolutionary rate) matrix. The number of bootstrap replicates
### is very small, 5 (to reduce running times). In reality it should be much more, e.g.,
### 100 or more if the computational budget allows.
BMbootstrap<-parametric.bootstrap(estimated.model=BMestim,phyltree=phyltree,
values.to.bootstrap=c("vX0","StS"),M.error=NULL,numboot=5)
### For the OUOU model, we want a CI for the evolutionary regression coefficient vector.
### The number of bootstrap replicates is very small, 5 (to reduce running times). 
### In reality it should be much more, e.g., 100 or more if the 
### computational budget allows.
OUOUbootstrap<-parametric.bootstrap(estimated.model=estimResultsOUOU,phyltree=phyltree,
values.to.bootstrap=c("evolutionary.regression"),regimes=regimes,root.regime="small",
M.error=NULL,predictors=c(3),kY=NULL,numboot=5,Atype=NULL,Syytype=NULL,diagA=NULL)
### For the OUBM model, we want CIs for the evolutionary and optimal regression 
### coefficient vectors. The number of bootstrap replicates is very small, 
### 5 (to reduce running times). In reality it should be much more, e.g.,
### 100 or more if the computational budget allows.
OUBMbootstrap<-parametric.bootstrap(estimated.model=OUBMestim,phyltree=phyltree,
values.to.bootstrap=c("evolutionary.regression","optimal.regression"),
regimes=regimes,root.regime="small",M.error=NULL,predictors=c(3),kY=2,
numboot=5,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")
}
}
