\name{mvslouchModel}
\alias{mvslouchModel}
\title{
Estimate parameters under a (multivariate) OUBM model of evolution
}
\description{
The \code{mvslouchModel} function uses maximum likelihood to fit parameters of a multivariate 
OUBM model evolving on the phylogeny.
}
\usage{
mvslouchModel(phyltree, data, kY, regimes = NULL, regimes.times = NULL, 
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
A data frame with the rows corresponding to the species while the columns correspond 
to the traits. The rows can be named by species, if not then the order
of the species has to be the same as the order in which the species are on the phylogeny.
}
  \item{kY}{
Number of "Y" (response) variables.
}
  \item{regimes}{
  A vector or list of regimes. If vector then each entry corresponds to the branch
  preceeding the respective node. If list then each list entry corresponds to a node
  and is a vector for regimes on that lineage. If \code{NULL} then a constant regime is 
  assumed on the whole tree.
}
  \item{regimes.times}{
A list of vectors for each tree node, it starts with 0 and ends with the current time 
of the species. In between are the times where the regimes (niches) changed. If \code{NULL} 
then each branch is considered to be a regime.
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
 An optional measurement error covariance matrix. The program tries to recognizes the 
 structure of matrix passed and accepts the following possibilities :
\itemize{
\item{}{a single number that will be on the diagonal of the covariance matrix}
\item{}{a m element vector with each value corresponding to a variable and the covariance matrix 
will have that vector repeated on its diagonal,}
\item{}{a nxm element vector a diagonal matrix with this vector on the diagonal,}
\item{}{a m x m ((number of variables) x (number of variables)) matrix it is assumed that the 
measurement errors are independent between observations so the resulting covariance structure is 
block diagonal,}
\item{}{a list of length m (number of variables), each list element is the covariance structure for 
the appropriate variable, either a single number (each observations has same variance), 
vector (of length n for each observation), or full matrix,}
\item{}{matrix of size mn x mn (m - number of variables, n - number of observations) 
the measurement error covaraince provided as is,}
\item{}{NULL no measurement error}
}
}


  \item{Atype}{
What class does the A matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"},
\code{"Symmetric"}, \code{"SymmetricPositiveDefinite"},
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
  An optional logical variable indicating whether to calculate and return confidence intervals 
  on the estimated parameters.
}
  \item{diagA}{
  Whether the values on \code{A}'s diagonal are to be \code{"Positive"}, \code{"Negative"} or sign allowed to vary, \code{NULL}.  
}
}
\details{
This function estimates the parameters of the following multivariate SDE,
\deqn{\begin{array}{rclccl}
dY(t) & = & -A(Y(t)-(\Psi(t)- A^{-1}BX(t)))dt + \Sigma_{yy} dB(t) & Y(0) & = & Y_{0}, \\ dX(t) & = & \Sigma_{xx} dB(t) & X(0) & = & X_{0}
\end{array}
}{
dY(t) = -A(Y(t)-(\Psi(t)- solve(A)BX(t)))dt + Syy dB(t) Y(0)=Y0, 

dX(t) = Sxx dB(t), X(0)=X0
}
on a phylogenetic tree. It uses
a numerical optimization over A (parametrized by its eigenvalues and eigenvectors or 
its QR decomposition)
and S (parametrized by its values) and conditional on A and S estimates the values of 
Psi corresponding
to the different regimes by a GLS estimate. Y(0) is assumed to be equal to - solve(A)BX(0) 
plus the root value of Psi.

The found point is described by a list containing four fields. 
The first field HeuristicSearchPointFinalFind is the parametrization of the model 
parameters at the considered
point with the value of the log--likelihood.
The field ParamsInModel is the point estimate of the parameters of the SDE. The value of mPsi0
can be currently ignored.
The field ParamSummary are different composite (evaluated at time 1) and summary statitics,
The field phylhalflife are the eigenvalues, eigenvectors and phylogenetic half lives
associated with the A matrix of, expmtA is
\eqn{exp(-A*1)}, optimal regression is the -solve(A)B matrix,
mPsi.rotated is each of the regime effects multiplied
by \eqn{1-exp(-A*1}, mPsi0.rotated is the same for the mPsi0 vector and can be ignored,
cov.matrix is the trait vector covariance matrix at time 1, 
corr.matrix is the trait vector correlation matrix at time 1, 
conditional.cov.matrix is the conditional covariance matrix of the OU type variables 
on the Brownian motion type at time \eqn{1}, i.e. Cov[Y|X](1),
conditional.corr.matrix is the conditional correlation matrix of the OU type variables 
on the Brownian motion type at time \eqn{1}, i.e. Corr[Y|X](1),
stationary.cov.matrix is the limit of the conditional.cov.matrix, 
stationary.corr.matrix is the limit of the conditional.corr.matrix, 
optima.cov.matrix
is the covariance matrix of the optimal process at time 1 equalling 1* solve(A)BSxx(Sxx^T)(B^T)solve(A)^T,
optima.corr.matrix
is the correlation matrix of the optimal process at time 1,
cov.with.optima is the covariance matrix between the optimal process and the Y type variables process,
corr.with.optima is the correlation matrix between the optimal process and the Y type variables process,
evolutionary.regression is the regression coefficient of E[Y|X](1)
StS the infinitesmal covariance matrix,
LogLik the log--likelihood, dof the degrees of freedom, m2loglik is -2log--likelihood,
aic is the Akaike information criterion, aic.c is the Akaike information criterion corrected for small
sample size, sic is the Schwarz information criterion, bic is the Bayesian information criterion
(which is the same as the Schwarz information criterion) and RSS is the residual sum of squares.
The last field LogLik is the log--likelihood at the point.
}
\value{
\item{FinalFound}{The point where the search procedure stopped. See Details for the description.}
\item{MaxLikFound}{The point with the highest likelihood found by the search procedure, 
if it is the same as the final point then this field equals "Same as final found". }
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

Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and 
  biological errors in phylogenetic comparative studies.
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
\note{
The slouch package is a recommended alternative if one has only a single response (Y) trait.
}

\section{Warning}{The estimation can take a long time and should be repeated a couple of times 
so that it is run from different starting positions.
The function can produce (a lot of) warnings and errors during the search procedure, 
this is nothing to worry about.
}

\seealso{
\code{slouch::model.fit}, \code{\link{SummarizeMVSLOUCH}}, \code{\link{simulMVSLOUCHProcPhylTree}}
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

### And finally try to recover the parameters of the mvOUBM model.
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive")
}
}
\keyword{models}

