## complex Hermitian matrices ("chm"); setClass("complex_herm_matrix") in 'aaa_allclasses.R'

`complex_herm_matrix` <- function(M){new("complex_herm_matrix",x=cbind(M))}  # this is the only place new("real_symmetric_matrix",...) is called
`is.complex_herm_matrix` <- function(x){inherits(x,"complex_herm_matrix")}

`r_to_n_chm` <- function(r){sqrt(r)}
`n_to_r_chm` <- function(n){n^2}

`is_ok_chm` <- function(r){ # 'r' = number of rows in [rowwise] matrix
    jj <- sqrt(r)
    if(jj == round(jj)){
        return(jj) # size of nxn complex hermitian matrix
    } else {
        stop("not correct")
    }
}

`valid_chm` <- function(object){
    x <- object@x
    if(!is.numeric(x)){
        return("not numeric")
    } else if(!is.matrix(x)){
        return("not a matrix")
    } else if(is_ok_chm(nrow(x)) < 0){
        return("must have appropriate size")
    } else {
        return(TRUE)
    }
}

setValidity("complex_herm_matrix", valid_chm)

`as.complex_herm_matrix` <- function(x,d,single=FALSE){  # single modelled on as.onion()
    if(is.complex_herm_matrix(x)){
        return(x)
    } else if(is.matrix(x)){
        return(complex_herm_matrix(x))
    } else if(is.vector(x)){
        if(single){
            return(complex_herm_matrix(x))
        } else {
            return(numeric_to_complex_herm_matrix(x,d)) # problem! we do not know how big it is
        }
    } else {
        stop("not recognised")
    }
}

`numeric_to_complex_herm_matrix` <- function(x,d){stop("no unique coercion")}

`rchm` <- function(n=3,d=5){complex_herm_matrix(matrix(round(rnorm(n*(d*d)),2),ncol=n))}
`chm_id` <- function(n,d){as.complex_herm_matrix(kronecker(chm1_to_vec(diag(nrow=d)),t(rep(1,n))))}

`vec_to_chm1` <- function(x){
   r <- length(x)
   n <- sqrt(r)
   stopifnot(n==round(n))
   out <- matrix(0i,n,n)
   out[upper.tri(out,FALSE)] <- x[(n+1):(n*(n+1)/2)] + 1i*x[(n*(n+1)/2+1):(n^2)]
   out <- out + ht(out)
   diag(out) <- x[seq_len(n)]
   return(out)  # complex hermitian matrix
}

`chm1_to_vec` <- function(M){
    c(
        Re(diag(M)),
        Re(M[upper.tri(M,FALSE)]),
        Im(M[upper.tri(M,FALSE)])
    )
}

`vec_chmprod_vec` <- function(x,y){
    x <- vec_to_chm1(x)
    y <- vec_to_chm1(y)
    chm1_to_vec((cprod(x,y)+cprod(y,x))/2)
}

setMethod("as.1matrix","complex_herm_matrix",function(x,drop=TRUE){
    out <- lapply(seq_along(x), function(i){x[i,drop=TRUE]})
    if((length(x)==1) & drop){out <- out[[1]]}
    return(out)
} )

`chm_prod_chm`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- vec_chmprod_vec(jj[[1]][,i],jj[[2]][,i])
    }
    return(as.jordan(out,e1))
}

`chm_inverse` <- function(e1){
    out <- as.matrix(e1)
    for(i in seq_len(ncol(out))){
        out[,i] <- chm1_to_vec(solve(e1[i,drop=TRUE])) # the meat
    }
    return(as.jordan(out,e1))
}

`chm_power_numeric` <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- chm1_to_vec(mymatrixpower(vec_to_chm1(jj[[1]][,i]),jj[[2]][i])) # the meat
    }
    return(as.jordan(out,e1))
}


`chm_arith_chm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_jordan(e1, e2),
         "-" = jordan_plus_jordan(e1,jordan_negative(e2)),
         "*" = chm_prod_chm(e1, e2),
         "/" = chm_prod_chm(e1, chm_inverse(e2)),
         "^" = stop("chm^chm not defined"),
         stop(paste("binary operator \"", .Generic, "\" not defined for chm"))
         )
}

`chm_arith_numeric` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e1, e2),  
         "-" = jordan_plus_numeric(e1,-e2),  
         "*" = jordan_prod_numeric(e1, e2),
         "/" = jordan_prod_numeric(e1, 1/e2),
         "^" = chm_power_numeric(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for chm"))
         )
}

`numeric_arith_chm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e2, e1),  
         "-" = jordan_plus_numeric(-e2,e1),  
         "*" = jordan_prod_numeric(e2, e1),
         "/" = jordan_prod_numeric(chm_inverse(e2),e1),
         "^" = jordan_power_jordan(e2, e1),
         stop(paste("binary operator \"", .Generic, "\" not defined for chm"))
         )
}

setMethod("Arith",signature(e1 = "complex_herm_matrix", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = jordan_negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on chm objects"))
                   )
          } )

setMethod("Arith",signature(e1="complex_herm_matrix",e2="complex_herm_matrix"),    chm_arith_chm    )
setMethod("Arith",signature(e1="complex_herm_matrix",e2="numeric"            ),    chm_arith_numeric)
setMethod("Arith",signature(e1="numeric"            ,e2="complex_herm_matrix"),numeric_arith_chm    )

setMethod("[",signature(x="complex_herm_matrix",i="index",j="missing",drop="logical"),
          function(x,i,j,drop){
              out <- as.matrix(x)[,i,drop=FALSE]
              if(drop){
                  if(ncol(out)==1){
                      return(vec_to_chm1(c(out)))
                  } else {
                      stop("for >1 element, use as.list()")
                  } 
              } else {
                  return(as.jordan(out,x))
              }
          } )
              
setReplaceMethod("[",signature(x="complex_herm_matrix",i="index",j="missing",value="complex_herm_matrix"),  # matches rsm equivalent
                 function(x,i,j,value){
                   out <- as.matrix(x)
                   out[,i] <- as.matrix(value)  # the meat
                   return(as.jordan(out,x))
                 } )

setReplaceMethod("[",signature(x="complex_herm_matrix",i="index",j="ANY",value="ANY"),function(x,i,j,value){stop("second argument redundant")})  # matches rsm equivalent
