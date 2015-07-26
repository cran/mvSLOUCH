\name{SummarizeBM}
\alias{SummarizeBM}
\title{Summarize parameters estimated under a Brownian motion model}
\description{Compiles a summary (appropriate moments, conditional moments, information criteria)
of parameters of a Brownian motion model at a given time point.
}
\usage{
SummarizeBM(phyltree, data, modelParams, t = c(1), dof = NULL, M.error = NULL, 
predictors = NULL, calcCI = FALSE)
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
A list of model parameters, as returned in \code{ParamsInModel} part of 
\cr \code{BrownianMotionModel}'s output.
}
  \item{t}{A vector of time points at which the summary is to be calculated. This allows for
  one to study (and plot) the (conditional) mean and covariance as functions of time.
}
  \item{dof}{Number of unknown parametees in the model, can be extracted from the output of 
  \code{BrownianMotionModel()}. If not provided all parameters are assumed unknown.
}
  \item{M.error}{
 An optional measurement error covariance matrix. The program tries to recognizes the structure of 
 matrix passed and accepts the following possibilities :
\itemize{
\item{}{a single number that will be on the diagonal of the covariance matrix}
\item{}{a m element vector with each value corresponding to a variable and the covariance matrix will 
have that vector repeated on its diagonal,}
\item{}{a nxm element vector a diagonal matrix with this vector on the diagonal,}
\item{}{a m x m ((number of variables) x (number of variables)) matrix it is assumed that the 
measurement errors are independent between observations so the resulting covariance structure is 
block diagonal,}
\item{}{a list of length m (number of variables), each list element is the covariance structure for 
the appropriate variable, either a single number (each observations has same variance), vector 
(of length n for each observation), or full matrix,}
\item{}{matrix of size mn x mn (m - number of variables, n - number of observations) the measurement 
error covaraince provided as is,}
\item{}{NULL no measurement error}
}
}


  \item{predictors}{
A vector giving the numbers of the columns from
\code{data} which are to be considered predictor ones, \emph{i.e.} conditioned
on in the program output.
}
  \item{calcCI}{
 An optional logical variable indicating whether to calculate and return confidence intervals on the 
 estimated parameters.
}
}

\value{A list for each provided time point. See the help of \code{\link{BrownianMotionModel}} for 
what the summary at each time point is.}

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
  Interpreting the evolutionary regression: the interplay between observational and biological 
  errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}
\author{Krzysztof Bartoszek}

\section{Warning}{Calculating the confidence intervals can take a very long time.
Warnings and errors (even of lot of them) can be produced during the confidence interval calculation, 
this is nothing to worry about.
}

\seealso{
\code{\link{BrownianMotionModel}}, \code{\link{simulBMProcPhylTree}}
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

### Define Brownian motion parameters to be able to simulate data 
### under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data and remove the values corresponding to the internal nodes.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### Recover the parameters of the Brownian motion.
BMestim<-BrownianMotionModel(phyltree,BMdata)

### And summarize them.
BM.summary<-SummarizeBM(phyltree,BMdata,BMestim$ParamsInModel,t=c(1),
dof=BMestim$ParamSummary$dof,calcCI=FALSE)
### if one would want the confidence intervals then set calcCI=TRUE
}
}
\keyword{models}
