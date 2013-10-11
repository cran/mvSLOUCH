\name{simulBMProcPhylTree}
\alias{simulBMProcPhylTree}
\title{
Simulate data on a phylogeny under a (multivariate) Brownian motion model
}
\description{
Simulate data on a phylogeny under a (multivariate) Brownian motion model
}
\usage{
simulBMProcPhylTree(phyltree, X0, Sigma, dropInternal = TRUE)
}
\arguments{
  \item{phyltree}{The phylogeny in \code{ouch} format. All of the internal nodes have to be uniquely named. The tree can be obtained from e.g. a \code{nexus} file by the \code{read.nexus} function
from the \code{ape} package and converted into the \code{ouch} format by \code{ouch}'s
\code{ape2ouch} function. See the example of how to correct the internal nodes.}
  \item{X0}{The ancestral, root state.}
  \item{Sigma}{The diffusion matrix of the Brownian motion.}
  \item{dropInternal}{Logical whether the simulated values at the internal nodes be changed to NA or not.}
}

\value{
Returns a data.frame with each row corresponding to a tree node and each column to a trait. If dropInternal is TRUE
then the entries for the internal nodes are changed to NAs.
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
\code{\link{BrownianMotionModel}}, \code{\link{SummarizeBM}}
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

### Define Brownian motion parameters to be able to simulate data 
### under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data and remove the values corresponding to the internal nodes.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
}	
}
\keyword{models}
