% --- Source file: fitch.mvsl.R ---
\name{fitch.mvsl}
\alias{fitch.mvsl}
\title{Unordered Fitch parsimony reconstruction of discrete character states}
\description{
Implements an unordered Fitch parsimony reconstruction of discrete niche     variables for use in the OU models where optima are modeled on discrete, categorical niche encodings. Allows for delayed and accelerated transformations to deal with ambiguities. 
Function was originally the \code{fitch} function from the \pkg{slouch} package.
}
\usage{
fitch.mvsl(phyltree, niche, deltran = FALSE, acctran = FALSE, root = NULL)
}
\arguments{
  \item{phyltree}{The phylogenetic tree in \pkg{ouch} format.}
  \item{niche}{The specific niche variable in the \pkg{slouch} data.frame to be 	reconstructed, entered as data.frame\code{$}niche }
  \item{deltran}{Implements a delayed transformation algorithm in order to 	deal with ambiguous nodes}
  \item{acctran}{Implements an accelerated transformation algorithm to deal 	with ambiguous nodes}
  \item{root}{An optional argument allowing the user to define a character 	state for the root (useful if the root node is ambiguously 		reconstructed)}
}
\details{
The \code{fitch.mvsl} function is meant to be interactive, where the user acts on the advice given in the 
returned messages whilst attempting to reconstruct ancestral states. 
If the root node is ambiguous after an initial reconstruction 
(a message will be printed to the screen if this is the case), 
this needs to be set by the user using the \code{root = "state"} argument in the function call. 
Any remaining ambiguous nodes can then be dealt with by specifying \code{deltran} or \code{acctran ="TRUE"} 
in the function call 
}
\value{
The \code{fitch.mvsl} function returns a vector of reconstructed character states.
}
\references{ 
Fitch, M.W. (1971) Defining the course of Evolution: Minimum change for a specific tree topology. \emph{Systematic Zoology} \bold{20}:406--416.

Swofford, D. L. and W.P. Maddison (1987) Reconstructing ancestral character states under Wagner parsimony. \emph{Mathematical Biosciences} \bold{87}: 199--229. 
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{
\url{slouch::fitch}, \url{slouch::slouchtree.plot}, \url{slouch::model.fit}, \url{slouch::ouch2slouch}
}

\examples{
phyltree<-ape2ouch(rtree(5))
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
1:(phyltree@nnodes-phyltree@nterm))

regimes<-c("A","B","B","C","C")
regimesFitch<-fitch.mvsl(phyltree,regimes,root=1,deltran=TRUE)
}
\keyword{models}
