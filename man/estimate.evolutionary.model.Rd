\name{estimate.evolutionary.model}
\alias{estimate.evolutionary.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function to find best (out of BM, OU, OUOU, OUBM) fitting evolutionary model and 
estimate its parameters.
}
\description{
The \code{estimate.evolutionary.model} function calls the
\code{BrownianMotionModel},  \code{ouchModel} and \code{mvslouchModel}
functions with different classes of evolutionary model parameters.
It then compares the resulting estimates by the AICc and 
returns the best overall model.
}
\usage{
estimate.evolutionary.model(phyltree, dfdata, regimes = NULL, 
root.regime = NULL, M.error = NULL, repeats = 3, model.setups = NULL, 
predictors = NULL, kY = NULL, doPrint = FALSE, pESS=NULL,estimate.root.state=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phyltree}{
The phylogeny in \code{ouch} format. All of the internal nodes have to be uniquely named.
The tree can be obtained from e.g. a \code{nexus} file by the \code{read.nexus} function
from the \code{ape} package and converted into the \code{ouch} format by \code{ouch}'s
\code{ape2ouch} function. See the example of how to correct the internal nodes.
It is strongly advisable to rescale the tree so that the determinant of the 
speciation time matrix (\code{phyltree@branch.times}) is about 1. This
significantly stabilizes the estimation procedure.
}
  \item{dfdata}{
A data frame with the rows corresponding to the species while the columns correspond to the 
traits.
The rows can be named by species, if not then the order
of the species has to be the same as the order in which the species are on the phylogeny.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to the branch
preceding the respective node. If list then each list entry corresponds to a node
and is a vector for regimes on that lineage. If \code{NULL} then a constant regime is assumed on
the whole tree. 
}
  \item{root.regime}{
The regime at the root of the tree. If not given taken as first element of regimes.
}
  \item{M.error}{
An optional measurement error covariance matrix. The program tries to recognizes the structure 
of matrix passed and accepts the following possibilities :
\itemize{
\item{}{a single number that will be on the diagonal of the covariance matrix}
\item{}{a m element vector with each value corresponding to a variable and the covariance matrix 
will have that vector repeated on its diagonal,}
\item{}{a nxm element vector a diagonal matrix with this vector on the diagonal,}
\item{}{a m x m ((number of variables) x (number of variables)) matrix it is assumed that the 
measurement errors are independent between observations so the resulting covariance structure is 
block diagonal,}
\item{}{a list of length m (number of variables), each list element is the covariance structure 
for the appropriate variable, either a single number (each observations has same variance), 
vector (of length n for each observation), or full matrix,}
\item{}{matrix of size mn x mn (m - number of variables, n - number of observations) the 
measurement error covariance provided as is,}
\item{}{NULL no measurement error}
}
}
  \item{repeats}{
How many starting points for the numerical maximum likelihood procedure should be tried for 
each model setup.
}

\item{model.setups}{
What models to try when searching for the best evolutionary model. 
This field may remain \code{NULL}, in this situation the function
generates using 

\code{.generate.basic.model.setups()} a basic list of models.
Allowed values are 
\itemize{
\item{\code{"basic"} }{A basic list of models to try out is generated, defined using 

\code{.generate.basic.model.setups()}. This list should be usually enough.}
\item{\code{"fundamental"} }{A slightly extended list of models to try out is generated, 
defined using \code{.generate.fund.model.setups()}. Compared to \code{"basic"} a few more 
models are added.}
\item{\code{"extended"} }{An extension of the \code{"fundamental"} list of models to try out. 
Defined using \code{.generate.ext.model.setups()} which at the moment calls \cr
\code{generate.model.setups()}.}
\item{\code{"all"} }{All possible models are generated, using 
\code{.generate.all.model.setups()}. Running it will take an intolerable amount of time.}
}
Alternatively the use is also free to provide their own list of models in this variable. 
Each element of the list is a list with fields.
\itemize{
\item{\code{evolmodel} }{The evolutionary model, it may take one of the three values \code{"BM"} 
(Brownian motion model), \code{"ouch"} (OUOU model), \code{"mvslouch"} \cr (OUBM model).}
\item{\code{Atype} }{The class of the \code{A} matrix, ignored if \code{evolmodel} equals 
\code{"BM"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"}, \code{"SymmetricPositiveDefinite"},
\code{"Symmetric"}, \code{"DecomposablePositive"}, \code{"DecomposableNegative"}, \cr
\code{"DecomposableReal"}, \code{"Invertible"}.
}
\item{\code{Syytype} }{The class of the \code{A} matrix, ignored if \code{evolmodel} equals 
\code{"BM"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"},  \code{"Symmetric"}, \code{"Any"}. 
}
\item{\code{diagA} }{Should the diagonal of \code{A} be forced to be positive (\code{"Positive"}),
negative (\code{"Negative"}) or the sign free to vary (\code{NULL})}
}
A minimum example list is \code{list(list(evolmodel="BM"))}.
}
  \item{predictors}{
A vector giving the numbers of the columns from
\code{dfdata} which are to be considered predictor ones, \emph{i.e.} conditioned
on in the program output. If not provided then the "X" variables 
are treated as predictors, but this only for the OUBM models (for the others in this case 
none are treated as predictors).
}
  \item{kY}{
Number of "Y" (response) variables, for the OUBM models.
}
  \item{doPrint}{
Should the function print out information on what it is doing (TRUE) or keep silent 
(default FALSE).
}
  \item{pESS}{
Should the function also find the best model taking into account the phylogenetic effective sample size
and it so what method. If \code{NULL} then do not take this into account. Otherwise one
of \code{"reg"} ("regression" effective sample size that takes into account all of the 
correlations between species explicitely),
\code{"mean"} (mean value effective sample size \eqn{1^{T}R^{-1}1}, where \eqn{R} 
is the interspecies correlation matrix),
\code{"MI"} (mutual information effective sample size), 
\code{"mvreg"}(multivariate version of "regression" effective sample size when each species is described by a suite of traits) , 
\code{"mvMI"} (multivariate mutual information effective sample size when each species is described by a suite of traits) indicating
the way to calculate the pESS. The default is \code{"reg"}.
}
\item{estimate.root.state}{
    Should the root state be estimate \code{TRUE} (not recommended) or set at the optimum 
    \code{FALSE} (recommended). Root state estimation is usually unreliable hence
    if fossil measurements are available prediction based on them and the estimated
    model will probably be more accurate. If there is only one regime, then 
    estimation of the root state separately is impossible and will not be allowed.
}
}
\details{
If \code{model.setups} is left at the default value the 
function will take a long time to run, as it performs estimation for each model 
(\code{generate.model.setups} generates 90 setups) times the value in repeats. 
Therefore if the user has particular hypotheses in mind then it is advisable to
prepare their own list.
}
\value{
A list is returned that describes the results of the search. See the help for
\code{\link{BrownianMotionModel}}, \code{\link{ouchModel}} and
\code{\link{mvslouchModel}} for the description of the lower level entries.
The elements of this list are the following
\item{\code{BestModel} }{The resulting best model found. Included are the model parameters, 
a "first-glance" qualitative description of the model, the most important parameters of the 
process (half-lives and regressions in the case of OU models) and what to call to obtain 
standard errors. It takes a long time to obtain them so calculating them is not part of the 
standard procedure.}
\item{\code{BestModelESS} }{Only if \code{pESS} was \code{TRUE}. The resulting best model found
taking into account the phylogenetic essential sample size. Included are the model parameters, 
a "first-glance" qualitative description of the model, the most important parameters of the 
process (half-lives and regressions in the case of OU models) and what to call to obtain 
standard errors. It takes a long time to obtain them so calculating them is not part of the 
standard procedure.}
\item{\code{testedModels} }{A list of results for each tried model.}
\item{\code{model.setups} }{A list of models tried.}
\item{\code{repeats} }{How many starting points were tried per model.} 
}

\references{
Ane, C. (2008) 
  Analysis of comparative data with hierarchical autocorrelation. 
  Annals of Applied Statistics 2:1078-1102.

Bartoszek, K. (2015)
  Phylogenetic effective sample size.
  ArXiv e-prints 1507.07113.

Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
  
Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

Hansen, T.F. and Pienaar, J. and Orzack, S.H. (2008) 
  A comparative method for studying adaptation to randomly evolving environment. 
  Evolution 62:1965-1977.
}
\author{Krzysztof Bartoszek}

\note{
The slouch package is a recommended alternative if one has a OUBM models and 
only a single response (Y) trait. 
The ouch package considers an univariate OU  model and looking at it could be helpful.
}

\seealso{
\code{\link{BrownianMotionModel}}, \code{\link{SummarizeBM}}, \code{\link{simulBMProcPhylTree}}
\code{\link[ouch]{hansen}}, \code{\link{ouchModel}}, \code{\link{SummarizeOUCH}}, 
\code{\link{simulOUCHProcPhylTree}},
\code{slouch::model.fit}, \code{\link{mvslouchModel}}, \code{\link{SummarizeMVSLOUCH}}, \cr
\code{\link{simulMVSLOUCHProcPhylTree}}
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

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data and remove the values corresponding to the internal nodes.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[-(1:(phyltree@nnodes-phyltree@nterm)),]

### And finally try to recover the parameters of the OUOU model.
estimResults<-estimate.evolutionary.model(demotreeouch,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),kY=2,
doPrint=TRUE,pESS=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
