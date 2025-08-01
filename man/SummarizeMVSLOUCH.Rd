\name{SummarizeMVSLOUCH}
\alias{SummarizeMVSLOUCH}
\title{
Summarize parameters estimated under a multivariate OUBM motion model
}
\description{
Compiles a summary (appropriate moments, conditional moments, information criteria)
of parameters of a multivariate OUBM model at a given time point. The user is 
recommended to install the suggested package \pkg{PCMBaseCpp} which significantly 
speeds up the calculations (see Details).
}
\usage{
SummarizeMVSLOUCH(phyltree, mData, modelParams, regimes = NULL, 
regimes.times = NULL, t = c(1), dof = NULL, M.error = NULL, predictors = NULL, 
Atype = "Invertible", Syytype = "UpperTri", min_bl = 0.0003, maxiter = 50)
}
\arguments{
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}. The \code{root.edge} field is ignored.
}
  \item{mData}{
A matrix with the rows corresponding to the tip species while the columns correspond to the traits.
The rows should be named by species \cr (field \code{phyltree$tip.label}), if not, then a warning
is thrown and the order of the species is assumed to be the same as the order in which the 
species are in the phylogeny (i.e. correspond to the node indices \code{1:n}, 
where \eqn{n}{n} is the number of tips). The columns should be named by traits, otherwise
a warning is thrown and generic names are generated.
}
  \item{modelParams}{
A list of model parameters, as returned in \code{ParamsInModel} part of \cr \code{mvslouchModel}'s output.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of \code{phyltree}'s
branches, i.e. to each row of \code{phyltree$edge}.If list then each list entry corresponds to a tip node
and is a vector for regimes on that lineage. If \code{NULL}, then a constant regime is assumed on
the whole tree.
}
  \item{regimes.times}{
A list of vectors for each tree node, it starts with 0 and ends with the current time of the species.
In between are the times where the regimes (niches) changed. If \code{NULL}, then each branch is
considered to be constant a regime.
}
  \item{t}{
A vector of time points at which the summary is to be calculated. This allows for
  one to study (and plot) the (conditional) mean and covariance as functions of time.
The function additionally returns the parameter summary at the tree's height.
}
  \item{dof}{
  Number of unknown parameters in the model, can be extracted from the output of
    \code{mvslouchModel()}. If not provided all parameters are assumed unknown.
}

  \item{M.error}{
An optional measurement error covariance structure. The measurement errors between species are assumed 
independent. The program tries to recognize the structure of the passed matrix and accepts the following 
possibilities :
\itemize{
\item{a single number that is a common measurement error for all tips and species,}
\item{a m element vector with each value corresponding to a variable, measurement errors are independent
between variables and each species is assumed to have the same measurement errors,}
\item{a m x m ((number of variables) x (number of variables)) matrix, all species will
have the same measurement error,}
\item{a list of length n (number of species), each list element is the covariance structure 
for the appropriate (numbering according to tree) species, either a single number 
(each variable has same variance), vector (of length m for each variable), or m x m matrix, 
the order of the list has to correspond to the order of the nodes in the \code{phyltree} object,}
\item{NULL no measurement error.}
}
From version \code{2.0.0} of \pkg{mvSLOUCH} it is impossible to pass a single joint measurement error matrix 
for all the species and traits. 
}
  \item{predictors}{
  A vector giving the numbers of the columns from
  \code{data} which are to be considered predictor ones, \emph{i.e.} conditioned
  on in the program output. If not provided then columns \code{(kY+1):ncol(mData)}, i.e. the 
"BM" ones, are treated as predictors.
}
  \item{Atype}{
What class does the A matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"}, 
\code{Symmetric}, \code{SymmetricPositiveDefinite}, \cr
\code{"DecomposablePositive"}, 
\code{"DecomposableNegative"}, \cr
\code{"DecomposableReal"}, \code{"Invertible"},
\code{"TwoByTwo"}, \code{"Any"}
}
  \item{Syytype}{
What class does the Syy matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"}, 
\code{"Symmetric"}, \code{"Any"}
}
  \item{min_bl}{Value to which \pkg{PCMBase}'s \code{PCMBase.Threshold.Skip.Singular} should be set.
  It indicates that branches of length shorter than \code{min_bl} should be skipped in likelihood
  calculations. Short branches can result in singular covariance matrices for the transition
  density along a branch. The user should adjust this value if a lot of warnings are raised by
  \pkg{PCMBase}  about singularities during the likelihood calculations.
  However, this does not concern tip branches-these cannot be skipped and hence should be long
  enough so that numerical issues are not raised.
  }
\item{maxiter}{The maximum number of iterations inside the GLS estimation procedure. In the first step 
  regime optima and \code{B} (and perhaps initial state) are estimated conditional on the other parameters 
  and current estimate   of \code{B}, then the estimate of \code{B} is update and the same phylogenetic GLS 
  is repeated.
  }

}

\details{
The likelihood calculations are done by the \pkg{PCMBase} package. However, there is a
C++ backend, \pkg{PCMBaseCpp}. If it is not available, then the likelihood is calculated
slower using pure R. However, with the calculations in C++ up to a 100-fold increase in
speed is possible (more realistically 10-20 times). The \pkg{PCMBaseCpp} package is available
from \url{https://github.com/venelin/PCMBaseCpp}.

The setting \code{Atype="Any"} means that one assumes the matrix \code{A} is eigendecomposable.
If \code{A} is defective, then the output will be erroneous.

From version \code{2.0.0} of \pkg{mvSLOUCH} the data has to be passed as a matrix.
To underline this the data parameter's name has been changed to \code{mData}. 

From version \code{2.0.0} of \pkg{mvSLOUCH} the parameter \code{calcCI} has been removed. 
The package now offers the possibility of bootstrap confidence intervals, see
function \code{parametric.bootstrap}.
}
\value{A list for each provided time point. See the help of \code{\link{mvslouchModel}} for what the 
summary at each time point is.}

\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.

Bartoszek, K. and Tredgett Clarke J. and Fuentes-Gonzalez, J. and Mitov, V. and Pienaar, J. 
  and Piwczynski, M. and Puchalka, R. and Spalik, K. and Voje, K. L. (2024) 
  Fast mvSLOUCH: Multivariate Ornstein-Uhlenbeck-based models of trait evolution on 
  large phylogenies, Methods in Ecology and Evolution 15(9):1507-1515.

Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

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
  Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and niche tracking. 
  The American Naturalist 174:204-220.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}
\author{Krzysztof Bartoszek}

\seealso{
\code{slouch::model.fit}, \code{\link{mvslouchModel}}, \code{\link{simulMVSLOUCHProcPhylTree}}, 
\code{\link{parametric.bootstrap}}
}
\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(3)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
regimes<-c("small","small","large","small")

### Define SDE parameters to be able to simulate data under the mvOUBM model.
OUBMparameters<-list(vY0=matrix(1,ncol=1,nrow=1),A=matrix(0.5,ncol=1,nrow=1),
B=matrix(c(2),ncol=1,nrow=1),mPsi=cbind("small"=1,"large"=-1),
Syy=matrix(2,ncol=1,nrow=1),vX0=matrix(0,ncol=1,nrow=1),Sxx=diag(1,1,1),
Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))

### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

## Here we do not do any recovery step
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMparameters,
regimes,t=c(1),dof=7,maxiter=2)

RNGversion(as.character(getRversion()))
\dontrun{ ##It takes too long to run this
## now less trivial simulation setup
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
regimes<-c("small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the mvOUBM model.
OUBMparameters<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,0),c(0,5)),
B=matrix(c(2,-2),ncol=1,nrow=2),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(0,1,1),Sxx=matrix(1,1,1),
Syx=matrix(0,ncol=1,nrow=2),Sxy=matrix(0,ncol=2,nrow=1))

### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,]

### Try to recover the parameters of the mvOUBM model.
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")

### Summarize them.
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMestim$FinalFound$ParamsInModel,
regimes,t=c(phyltree$tree_height),dof=OUBMestim$FinalFound$ParamSummary$dof,maxiter=50)

### And finally bootstrap with particular interest in the evolutionary and optimal
### regressions
OUBMbootstrap<-parametric.bootstrap(estimated.model=OUBMestim,phyltree=phyltree,
values.to.bootstrap=c("evolutionary.regression","optimal.regression"),
regimes=regimes,root.regime="small",M.error=NULL,predictors=c(3),kY=2,
numboot=5,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")

}
}
\keyword{models}
\keyword{multivariate}
\keyword{htest}
