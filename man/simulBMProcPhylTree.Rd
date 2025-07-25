\name{simulBMProcPhylTree}
\alias{simulBMProcPhylTree}
\title{
Simulate data on a phylogeny under a (multivariate) Brownian motion model
}
\description{
Simulate data on a phylogeny under a (multivariate) Brownian motion model
}
\usage{
simulBMProcPhylTree(phyltree, X0, Sigma, dropInternal = TRUE, M.error=NULL, 
fullTrajectory=FALSE, jumpsetup=NULL, keep_tree = FALSE, step=NULL)
}
\arguments{
  \item{phyltree}{
  The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a 
  \code{nexus} file by the \code{read.nexus()} function from the \pkg{ape} package. 
  The "standard" \pkg{ape} node indexing is assumed: for a tree with \eqn{n}{n} tips, 
  the tips should have indices \code{1:n} and the root index \code{n+1}.
  The \code{root.edge} field is ignored.
  }
  \item{X0}{The ancestral, root state.}
  \item{Sigma}{The diffusion matrix of the Brownian motion.}
  \item{dropInternal}{Logical whether the simulated values at the internal nodes should be 
  dropped.}
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
    \item{fullTrajectory}{Should the full realization of the process or only node and tip values 
    be returned}
    \item{jumpsetup}{Either \code{NULL} or list describing the jump at speciation. In the second 
    case:
    \describe{
	\item{\code{jumptype}}{In what way does the jump take place. Possible values are \cr
	\code{"ForBoth"} the jump occurs at speciation and is common to both daughter lineages,
	\code{"RandomLineage"} the jump occurs just after speciation affecting exactly one
	daughter lineage, both descending branches have the same chance of being affected,
	\code{"JumpWithProb"} the jump occurs with probability \code{jumpprob}
	just after speciation independently on each daughter lineage independently.
	}
	\item{\code{jumpprob}}{A value in \eqn{[0,1]} indicating the probability of a 
	jump taking place, only matters if \code{jumptype} is \code{"JumpWithProb"} 
	 or \code{"JumpWithProb"}.}
	\item{\code{jumpdistrib}}{The distribution of the jump, currently only can take
	value \code{"Normal"}.}
	\item{\code{vMean}}{The expected value of the jump, a vector of appropriate
	length if the trait is multivariate.}
	\item{\code{mCov}}{The variance of the jump, a matrix of appropriate
	dimensions if the trait is multivariate.}
    }
    }
    \item{keep_tree}{Logical whether the used tree should be saved inside the output object.
    Useful for any future reference, but as the tree is enhanced for \pkg{mvSLOUCH}'s needs
    the resulting output object may be very large (it the number of tips is large).}
    \item{step}{The step size of the simulation.}
}

\value{
If \code{fullTrajectory} is \code{FALSE} then 
returns a matrix with each row corresponding to a tree node and each column to a trait. 
Otherwise returns a more complex object describing the full realization of the process on the tree.
If \code{dropInternal} is \code{TRUE}, then the entries for the internal nodes are changed to 
\code{NA}s. The ordering of the rows corresponds to the order of the nodes (their indices) in 
the \code{phylo} object. Hence, the first \code{n} rows will be the tip rows 
(by common \code{phylo} convention).
}
\references{
Bartoszek, K. (2014) 
  Quantifying the effects of anagenetic and cladogenetic evolution. 
  Mathematical Biosciences 254:42-57.

Bartoszek, K. (2016) 
  A Central Limit Theorem for punctuated equilibrium.
  arXiv:1602.05189.

Bartoszek, K. and Fuentes-Gonzalez, J. and Mitov, V. and Pienaar, J. and Piwczynski, M.
  and Puchalka, R. and Spalik, K. and Voje, K. L. (2023) 
  Model Selection Performance in Phylogenetic Comparative Methods Under Multivariate 
  Ornstein-Uhlenbeck Models of Trait Evolution, Systematic Biology 72(2):275-293.

Bartoszek, K. and Fuentes-Gonzalez, J. and Mitov, V. and Pienaar, J. and Piwczynski, M.
  and Puchalka, R. and Spalik, K. and Voje, K. L. (2024) 
  Analytical advances alleviate model misspecification in non-Brownian multivariate
  comparative methods, Evolution 78(3):389-400.

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

\seealso{
\code{\link{BrownianMotionModel}}, \code{\link{SummarizeBM}}
}
\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define Brownian motion parameters to be able to simulate data 
### under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data.
jumpobj<-list(jumptype="RandomLineage",jumpprob=0.5,jumpdistrib="Normal",
vMean=rep(0,3),mCov=diag(1,3,3))
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx,
jumpsetup=jumpobj)
RNGversion(as.character(getRversion()))

}
\keyword{models}
\keyword{multivariate}
\keyword{datagen}
