\name{SummarizeBM}
\alias{SummarizeBM}
\title{Summarize parameters estimated under a Brownian motion model}
\description{Compiles a summary (appropriate moments, conditional moments, information criteria)
of parameters of a Brownian motion model at a given time point. The user is recommended to install 
suggested package \pkg{PCMBaseCpp} which significantly speeds up the calculations (see Details).
}
\usage{
SummarizeBM(phyltree, mData, modelParams, t = c(1), dof = NULL, M.error = NULL, 
predictors = NULL, min_bl = 0.0003)
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
A list of model parameters, as returned in \code{ParamsInModel} part of 
\cr \code{BrownianMotionModel}'s output.
}
  \item{t}{A vector of time points at which the summary is to be calculated. This allows for
  one to study (and plot) the (conditional) mean and covariance as functions of time.
  The function additionally returns the parameter summary at the tree's height.  
}
  \item{dof}{Number of unknown parameters in the model, can be extracted from the output of 
  \code{BrownianMotionModel()}. If not provided all parameters are assumed unknown.
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
on in the program output. If not provided, then none will be treated as predictors.
}
  \item{min_bl}{Value to which \pkg{PCMBase}'s \code{PCMBase.Threshold.Skip.Singular} should be set.
  It indicates that branches of length shorter than \code{min_bl} should be skipped in likelihood
  calculations. Short branches can result in singular covariance matrices for the transition
  density along a branch. The user should adjust this value if a lot of warnings are raised by
  \pkg{PCMBase}  about singularities during the likelihood calculations.
  However, this does not concern tip branches-these cannot be skipped and hence should be long
  enough so that numerical issues are not raised.
  }
}

\details{
The likelihood calculations are done by the \pkg{PCMBase} package. However, there is a
C++ backend, \pkg{PCMBaseCpp}. If it is not available, then the likelihood is calculated
slower using pure R. However, with the calculations in C++ up to a 100-fold increase in
speed is possible (more realistically 10-20 times). The \pkg{PCMBaseCpp} package is available
from \url{https://github.com/venelin/PCMBaseCpp}.

The \code{phyltree_paths()} function enhances the tree for usage by \code{mvSLOUCH}. 
Hence, to save time, it is advisable to first do \code{phyltree<-mvSLOUCH::phyltree_paths(phyltree)} 
and only then use it with \code{BrownianMotionModel()}.

From version \code{2.0.0} of \pkg{mvSLOUCH} the data has to be passed as a matrix.
To underline this the data parameter's name has been changed to \code{mData}. 

From version \code{2.0.0} of \pkg{mvSLOUCH} the parameter \code{calcCI} has been removed. 
The package now offers the possibility of bootstrap confidence intervals, see
function \code{parametric.bootstrap}.
}
\value{A list for each provided time point. See the help of \code{\link{BrownianMotionModel}} for 
what the summary at each time point is.}

\references{
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
\code{\link{BrownianMotionModel}}, \code{\link{simulBMProcPhylTree}}, \code{\link{parametric.bootstrap}}
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
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[phyltree$tip.label,,drop=FALSE]

### Recover the parameters of the Brownian motion.
BMestim<-BrownianMotionModel(phyltree,BMdata)

### Summarize them.
BM.summary<-SummarizeBM(phyltree,BMdata,BMestim$ParamsInModel,t=c(1),
dof=BMestim$ParamSummary$dof)
RNGversion(as.character(getRversion()))
#\dontrun
{ ##It takes too long to run this
### Now obtain bootstrap confidence intervals for some parameters.
BMbootstrap<-parametric.bootstrap(estimated.model=BMestim,phyltree=phyltree,
values.to.bootstrap=c("vX0","StS"),,M.error=NULL,numboot=5)
}
}
\keyword{models}
\keyword{multivariate}
\keyword{htest}
