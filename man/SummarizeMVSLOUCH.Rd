\name{SummarizeMVSLOUCH}
\alias{SummarizeMVSLOUCH}
\title{
Summarize parameters estimated under a multivariate OUBM motion model
}
\description{
Compiles a summary (appropriate moments, conditional moments, information criteria)
of parameters of a multivariate OUBM model at a given time point.
}
\usage{
SummarizeMVSLOUCH(phyltree, data, modelParams, regimes = NULL, 
regimes.times = NULL, t = c(1), dof = NULL, M.error = NULL, predictors = NULL, 
Atype = "Invertible", Syytype = "UpperTri", calcCI = FALSE)
}
\arguments{
  \item{phyltree}{
The phylogeny in \code{ouch} format. All of the internal nodes have to be uniquely named.
The tree can be obtained from e.g. a \code{nexus} file by the \code{read.nexus} function
from the \code{ape} package and converted into the \code{ouch} format by \code{ouch}'s
\code{ape2ouch} function. See the example of how to correct the internal nodes.
}
  \item{data}{
A data frame with the rows corresponding to the species while the columns correspond to the traits.
The rows can be named by species, if not then the order
of the species has to be the same as the order in which the species are on the phylogeny.
}
  \item{modelParams}{
A list of model parameters, as returned in \code{ParamsInModel} part of \cr \code{mvslouchModel}'s output.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to the branch
preceeding the respective node. If list then each list entry corresponds to a node
and is a vector for regimes on that lineage. If \code{NULL} then a constant regime is assumed on
the whole tree.
}
  \item{regimes.times}{
A list of vectors for each tree node, it starts with 0 and ends with the current time of the species.
In between are the times where the regimes (niches) changed. If \code{NULL} then each branch is
considered to be a regime.
}
  \item{t}{
A vector of time points at which the summary is to be calculated. This allows for
  one to study (and plot) the (conditional) mean and covariance as functions of time.
}
  \item{dof}{
  Number of unknown parametees in the model, can be extracted from the output of
    \code{mvslouchModel()}. If not provided all parameters are assumed unknown.
}

  \item{M.error}{
 An optional measurement error covariance matrix. The program tries to recognizes the structure of matrix passed and accepts the following possibilities :
\itemize{
\item{}{a single number that will be on the diagonal of the covariance matrix}
\item{}{a m element vector with each value corresponding to a variable and the covariance matrix will have that vector repeated on its diagonal,}
\item{}{a nxm element vector a diagonal matrix with this vector on the diagonal,}
\item{}{a m x m ((number of variables) x (number of variables)) matrix it is assumed that the measurement errors are independent between observations so the resulting covariance structure is block diagonal,}
\item{}{a list of length m (number of variables), each list element is the covariance structure for the appropriate variable, either a single number (each observations has same variance), vector (of length n for each observation), or full matrix,}
\item{}{matrix of size mn x mn (m - number of variables, n - number of observations) the measurement error covaraince provided as is,}
\item{}{NULL no measurement error}
}
}

  \item{predictors}{
  A vector giving the numbers of the columns from
  \code{data} which are to be considered predictor ones, \emph{i.e.} conditioned
  on in the program output.
}
  \item{Atype}{
What class does the A matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"}, 
\code{Symmetric}, \code{SymmetricPositiveDefinite}, 
\code{"DecomposablePositive"}, \cr
\code{"DecomposableNegative"}, 
\code{"DecomposableReal"}, \code{"Invertible"},
\code{"TwoByTwo"}
}
  \item{Syytype}{
What class does the Syy matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"}, 
\code{Symmetric}, \code{Any}
}
  \item{calcCI}{
   An optional logical variable indicating whether to calculate and return confidence intervals on the estimated
parameters.
}
}

\details{If calcCI is set to TRUE the function returns the confidence intervals, the log-likelihood
surface for the eigenvalues conditional on the other parameters and regression confidence intervals.
See Bartoszek et. al. for the mathematical details of these confidence intervals.
}
\value{A list for each provided time point. See the help of \code{\link{mvslouchModel}} for what the summary at each time point is.}

\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.

Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

Hansen, T.F. (1997) 
  Stabilizing selection and the comparative analysis of adaptation. 
  Evolution 51:1341-1351.

Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
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

\section{Warning}{Calculating the confidence intervals can take a very long time.
Warnings and errors (even of lot of them) can be produced during the confidence interval calculation, this is nothing to worry about.
}

\seealso{
\code{slouch::model.fit}, \code{\link{mvslouchModel}}, \code{\link{simulMVSLOUCHProcPhylTree}}
}
\examples{
\dontrun{ ##It takes too long to run this
### We will first simulate a small phylogenetic tree using functions from ape and ouch.
### For simulating the tree one could also use alternative functions, eg. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape2ouch(rtree(5))

### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
1:(phyltree@nnodes-phyltree@nterm))

### Define a vector of regimes.
regimes<-c("small","small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the mvOUBM model.
OUBMparameters<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,0),c(0,5)),
B=matrix(c(2,-2),ncol=1,nrow=2),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(0,1,1),Sxx=matrix(1,1,1),
Syx=matrix(0,ncol=1,nrow=2),Sxy=matrix(0,ncol=2,nrow=1))

### Now simulate the data and remove the values corresponding to the internal nodes.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### Recover the parameters of the mvOUBM model.
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")

### And summarize them.
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMestim$FinalFound$ParamsInModel,
regimes,t=c(1),dof=OUBMestim$FinalFound$ParamSummary$dof,calcCI=FALSE)
### if one would want the confidence intervals then set calcCI=TRUE
}
}
\keyword{models}
