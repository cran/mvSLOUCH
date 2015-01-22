\name{BrownianMotionModel}
\alias{BrownianMotionModel}
\title{Estimate parameters under a Brownian motion model of evolution}
\description{
The \code{BrownianMotionModel} function uses maximum likelihood to fit parameters of a Brownian
motion model evolving on the phylogeny. If there are no missing values
or measurement error then it just calls the \code{brown()} function from the \code{ouch} package.
}
\usage{
BrownianMotionModel(phyltree, data, predictors = NULL, M.error = NULL, 
calcCI = FALSE)
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
  \item{predictors}{
A vector giving the numbers of the columns from
\code{data} which are to be considered predictor ones, \emph{i.e.} conditioned
on in the program output. 
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
  \item{calcCI}{
  An optional logical variable indicating whether to calculate and return confidence intervals on the estimated parameters.
}
}
\details{
This function estimates the parameters of a multivariate Brownian motion model defined by the SDE,
\deqn{dX(t) = \Sigma dB(t), X(0)=X_{0}}{dX(t) = S dW(t), X(0)=X0} evolving on a phylogenetic tree. If there is no measurement error it uses the \code{brown} function
from the ouch package, otherwise it does a numerical optimization.
}
\value{
    \item{ParamsInModel}{
    A list with estimated model parameters. The elements are vX0 : the ancestral trate, and Sxx where t*Sxx%*%t(Sxx) is the Brownian motion's covariance matrix at time t.
}
    \item{ParamSummary}{
    A list with summary statistics with elements,
    StS the infinitesmal covariance matrix Sxx%*%t(Sxx),
    LogLik the log--likelihood, dof the degrees of freedom, m2loglik is $-2$log--likelihood,
    aic is the Akaike information criterion, aic.c is the Akaike information criterion corrected for small
    sample size, sic is the Schwarz information criterion, bic is the Bayesian information criterion
    (which is the same as the Schwarz information criterion) and RSS is the residual sum of squares.
}
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
  
Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}

\author{Krzysztof Bartoszek}

\seealso{
\code{\link{SummarizeBM}}, \code{\link{simulBMProcPhylTree}}
}
\examples{
### We will first simulate a small phylogenetic tree  using functions from ape and ouch.
### For simulating the tree one could also use alternative functions, eg. sim.bd.taxa 
### from the TreeSim package
\dontrun{ ##It takes too long to run this
phyltree<-ape2ouch(rtree(5))

### Correct the names of the internal node labels
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
1:(phyltree@nnodes-phyltree@nterm))

### Define Brownian motion parameters to be able to simulate data under 
### the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data and remove the values corresponding to the internal nodes.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### And finally try to recover the parameters of the Brownian motion.
BMestim<-BrownianMotionModel(phyltree,BMdata)
}
}
\keyword{models}
