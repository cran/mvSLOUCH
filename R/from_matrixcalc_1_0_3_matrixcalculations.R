## This file is included in the mvSLOUCH R software package.

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## Krzysztof Bartoszek:
## This file is copied from matrixcalc_1.0.3 found in CRAN archive.
## The matrixcalc package was removed from CRAN and hence its functionality stopped being available via CRAN.
## The required files  
## (is.diagonal.matrix.R, is.positive.definite.R, is.positive.semi.definite.R, is.singular.matrix.R, is.square.matrix.R, is.symmetric.matrix.R) 
## from matrixcalc were copied here with only minor changes 
## 	function names have .matrixcalc_ added in front of them to indicate their origin
##	is.square.matrix() was changed .matrixcalc_is.square.matrix()
## matrixcalc was originally on CRAN licensed as  GPL-2 | GPL-3 [expanded from: GPL (>= 2)]
## and this is of course retained
## =========================================================================================================


.matrixcalc_is.diagonal.matrix <- function( x, tol=1e-8 )
{
## called in evolmodelest.R
###
### this function returns TRUE if the off diagonal elements in absolute
### value are less than the given tolerance and FALSE otherwise
###
### arguments
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!.matrixcalc_is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    y <- x
    diag( y ) <- rep( 0, nrow( y ) )
    return( all( abs( y ) < tol ) )
}

.matrixcalc_is.positive.semi.definite <- function( x, tol=1e-8 )
{
## called in PhyloSDE.R
###
### this function determines if the given real symmetric matrix is positive semi definite
### parameters
### x = a square numeric matrix object
### tol = tolerance level for zero
###
    if ( !.matrixcalc_is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !.matrixcalc_is.symmetric.matrix( x ) )
        stop( "argument x is not a symmetric matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    eigenvalues <- eigen(x, only.values = TRUE)$values
    n <- nrow( x )
    for ( i in 1: n ) {
        if ( abs( eigenvalues[i] ) < tol ) {
            eigenvalues[i] <- 0
        }
    }    
    if ( any( eigenvalues < 0 ) ) {
        return( FALSE )
    }
    return( TRUE )
}

.matrixcalc_is.singular.matrix <- function( x, tol=1e-8 )
{
###
### this function returns TRUE if the matrix argument x
### is singular and FALSE otherwise
###
### argument
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!.matrixcalc_is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    det.x <- det( x )
    return( abs( det.x ) < tol )
}

.matrixcalc_is.symmetric.matrix <- function( x )
{
## called in PhyloSDE.R, estimBM.R, evolmodelest.R, getESS.R, loglik.R, sdecovariancephyl.R, sdemoments.R
###
### this function determines if the matrix is symmetric
###
### argument
### x = a numeric matrix object
###
    if ( !is.matrix( x ) ) {
        stop( "argument x is not a matrix" )
    }
    if ( !is.numeric( x ) ) {
        stop( "argument x is not a numeric matrix" )
    }    
    if ( !.matrixcalc_is.square.matrix( x ) )
        stop( "argument x is not a square numeric matrix" )
    return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
}

.matrixcalc_is.positive.definite <- function( x, tol=1e-8 )
{
## called in estimBM.R, evolmodelest.R, getESS.R, loglik.R, sdecovariancephyl.R, sdemoments.R
###
### this function determines if the given real symmetric matrix is positive definite
###
### parameters
### x = a square numeric matrix object
### tol = tolerance level for zero
###
    if ( !.matrixcalc_is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !.matrixcalc_is.symmetric.matrix( x ) )
        stop( "argument x is not a symmetric matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    eigenvalues <- eigen(x, only.values = TRUE)$values
    n <- nrow( x )
    for ( i in 1: n ) {
        if ( abs( eigenvalues[i] ) < tol ) {
            eigenvalues[i] <- 0
        }
    }    
    if ( any( eigenvalues <= 0 ) ) {
        return( FALSE )
    }
    return( TRUE )
}

.matrixcalc_is.square.matrix <- function( x )
{
###
### determines if the given matrix is a square matrix
###
### arguments
### x = a matrix object
###
    if ( !is.matrix( x ) )
        stop( "argument x is not a matrix" )
    return( nrow(x) == ncol(x) )
}
