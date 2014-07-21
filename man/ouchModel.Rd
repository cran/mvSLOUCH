\name{ouchModel}
\alias{ouchModel}
\title{
Estimate parameters under a (multivariate) OU model of evolution
}
\description{
The \code{ouchModel} function uses maximum likelihood to fit parameters of a multivariate OU  model evolving on the phylogeny.

}
\usage{
ouchModel(phyltree, data, regimes = NULL, regimes.times = NULL, 
root.regime = NULL, predictors = NULL, M.error = NULL, Atype = "Invertible", 
Syytype = "UpperTri", calcCI = FALSE, diagA = "Positive")
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
  \item{root.regime}{
The regime at the root of the tree. If not given taken as first element of regimes.
}
  \item{predictors}{
A vector giving the numbers of the columns from
  \code{data} which are to be considered predictor ones, \emph{i.e.} conditioned
      on in the program output. If not provided the "X" variables are treated as predictors.
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


  \item{Atype}{
What class does the A matrix in the multivariate OUOU model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"},
\code{Symmetric}, \code{SymmetricPositiveDefinite}, 
\code{"DecomposablePositive"}, \cr
\code{"DecomposableNegative"}, \code{"DecomposableReal"}, \code{"Invertible"},
\code{"TwoByTwo"}
}
  \item{Syytype}{
What class does the Syy matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"},
\code{Symmetric}, \code{Any}
}
  \item{calcCI}{
  An optional logical variable indicating whether to calculate and return confidence intervals on the estimated parameters.
}
\item{diagA}{
  Whether the values on \code{A}'s diagonal are to be \code{"Positive"}, \code{"Negative"} or sign allowed to vary, \code{NULL}.
}
  
}
\details{
This function estimates the parameters of the following multivariate SDE,
\deqn{dY(t) = -A(Y-\Psi(t))dt + \Sigma dW(t), Y(0)=Y_{0}}{dY(t) = -A(Y-Psi(t))dt + S dW(t), Y(0)=Y0} on a phylogenetic tree. It uses
a numerical optimization over A (parametrized by its eigenvalues and eigenvectors or its QR decomposition)
and S (parametrized by its values) and conditional on A and S estimates the values of Psi corresponding
to the different regimes by a GLS estimate. Y(0) is assumed to be equal to the root value of Psi.

The found point is described by a list containing four fields. 
The first field HeuristicSearchPointFinalFind is the parametrization of the model parameters at the considered
point with the value of the log--likelihood.
The field ParamsInModel is the point estimate of the parameters of the SDE. The value of mPsi0
can be currently ignored.
The field ParamSummary are different composite (evaluated at time 1) and summary statitics,
The field phylhalflife are the eigenvalues, eigenvectors and phylogenetic half lives
associated with the A matrix, expmtA is exp(-A*1),
mPsi.rotated is each of the regime effects multiplied
by (1-exp(-A*1), mPsi0.rotated is the same for the mPsi0 vector and can be ignored,
cov.matrix is the trait vector covariance matrix at time 1, 
corr.matrix is the trait vector correlation matrix at time 1,
trait.regression
is a list consisting of regression coefficents when taking each trait in turn and
calculating its coniditional expectation on all of the other trait, 
stationary.cov.matrix is the stationary covariance matrix of process if it exists
(i.e. the eigenvalues have positive real part),
stationary.corr.matrix is the stationary correlation matrix of process if it exists
(i.e. the eigenvalues have positive real part),
StS the infinitesmal covariance matrix S%*%t(S),
LogLik the log--likelihood, dof the degrees of freedom, m2loglik is $-2$log--likelihood,
aic is the Akaike information criterion, aic.c is the Akaike information criterion corrected for small
sample size, sic is the Schwarz information criterion, bic is the Bayesian information criterion
(which is the same as the Schwarz information criterion) and RSS is the residual sum of squares.
The last field LogLik is the log--likelihood at the point.
}
\value{
\item{FinalFound}{The point where the search procedure stopped. See Details for the description.}
\item{MaxLikFound}{The point with the highest likelihood found by the search procedure, if it is 
the same as the final point then this field equals "Same as final found". }
}
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

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.  
}
\author{Krzysztof Bartoszek}
\note{
The ouch package consideres a similar model and looking at it could be helpful.
}

\section{Warning}{The estimation can take a long time and should be repeated 
a couple of times so that it is run from different starting positions.
The function can produce (a lot of) warnings and errors during the search procedure, this is nothing to worry about.
}

\seealso{
\code{\link[ouch]{hansen}}, \code{\link{SummarizeOUCH}}, \code{\link{simulOUCHProcPhylTree}}
}

\examples{
\donttest{
### We will first simulate a small phylogenetic tree using functions from ape and ouch.
### For simulating the tree one could also use alternative functions, eg. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape2ouch(rtree(5))

### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
1:(phyltree@nnodes-phyltree@nterm))

### Define a vector of regimes.
regimes<-c("small","small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data and remove the values corresponding to the internal nodes.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### And finally try to recover the parameters of the OUOU model.
OUOUestim<-ouchModel(phyltree,OUOUdata,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")
}
}
\keyword{models}

