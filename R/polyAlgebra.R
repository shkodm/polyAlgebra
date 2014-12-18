#' Package for analitical calculations on multiple variable polynomials
#'
#' @docType package
#' @name polyAlgebra-package
#' @rdname polyAlgebra-package
#' @useDynLib polyAlgebra
NULL

#' Class that represents a polynomial function
#' 
#' @slot tab Array keeping all the information about the function (see below)
#' @exportClass pAlg
pAlg = setClass("pAlg", representation(tab="array"))
