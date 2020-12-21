## real symmetric matrices; setClass("real_symmetric_matrix") is in  aaa_allclasses.R

`real_symmetric_matrix` <- function(M){new("real_symmetric_matrix",x=cbind(M))}  # this is the only place new("real_symmetric_matrix",...) is called
`is.real_symmetric_matrix` <- function(x){inherits(x,"real_symmetric_matrix")}

`is_ok_rsm` <- function(r){ # 'r' = number of rows in [rowwise] matrix
    jj <- sqrt(1+8*r)
    if(jj == round(jj)){
        return((jj+1)/2) # size of nxn real matrix
    } else {
        stop("not correct")
    }
}

`valid_rsm` <- function(object){
    x <- object@x
    if(!is.numeric(x)){
        return("not numeric")
    } else if(!is.matrix(x)){
        return("not a matrix")
    } else if(is_ok_rsm(nrow(x)) < 0){
        return("must have appropriate size")
    } else {
        return(TRUE)
    }
}

setValidity("real_symmetric_matrix", valid_rsm)

`as.real_symmetric_matrix` <- function(x,d,single=FALSE){  # single modelled on as.onion()
    if(is.real_symmetric_matrix(x)){
        return(x)
    } else if(is.matrix(x)){
        return(real_symmetric_matrix(x))
    } else if(is.vector(x)){
        if(single){
            return(real_symmetric_matrix(x))
        } else {
            return(scalars_to_real_symmetric_matrix(x,d)) # problem! we do not know how big it is
        }
    } else {
        stop("not recognised")
    }
}

`rrsm` <- function(n=3){real_symmetric_matrix(matrix(round(rnorm(n*10),2),nrow=10))}

`vec_to_rsm1` <- function(x){
   r <- length(x)
   n <- (sqrt(1+8*r)-1)/2
   stopifnot(n==round(n))
   out <- matrix(0,n,n)
   out[upper.tri(out,TRUE)] <- x
   out <- out + t(out)
   diag(out) <- diag(out)/2
   return(out)
}    

`vec_rsmprod_vec` <- function(x,y){
    x <- vec_to_rsm1(x)
    y <- vec_to_rsm1(y)
    jj <- (cprod(x,y)+cprod(y,x))/2
    return(jj[upper.tri(jj,TRUE)])
}

setMethod("as.1matrix","real_symmetric_matrix",function(x,drop=TRUE){
    out <- lapply(seq_along(x), function(i){x[i,drop=TRUE]})
    if((length(x)==1) & drop){out <- out[[1]]}
    return(out)
} )

`rsm_prod_rsm`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- vec_rsmprod_vec(jj[[1]][,i],jj[[2]][,i])
    }
    return(as.jordan(out,e1))
}

`rsm_arith_rsm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_jordan(e1, e2),
         "-" = jordan_plus_jordan(e1,jordan_negative(e2)),
         "*" = rsm_prod_rsm(e1, e2),
         "/" = rsm_prod_rsm(e1, jordan_matrix_inverse(e2)), # fails
         "^" = stop("rsm^rsm not defined"),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

`rsm_arith_numeric` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e1, e2),  
         "-" = jordan_plus_numeric(e1,-e2),  
         "*" = jordan_prod_numeric(e1, e2),
         "/" = jordan_prod_numeric(e1, 1/e2),
         "^" = rsm_power_numeric(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

`numeric_arith_rsm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e2, e1),  
         "-" = jordan_plus_numeric(-e2,e1),  
         "*" = jordan_prod_numeric(e2, e1),
         "/" = jordan_prod_numeric(e2, 1/e1),
         "^" = jordan_power_jordan(e2, e1),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

setMethod("Arith",signature(e1 = "real_symmetric_matrix", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = jordan_negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on alberts"))
                   )
          } )

setMethod("Arith",signature(e1="real_symmetric_matrix" ,e2="real_symmetric_matrix"),    rsm_arith_rsm    )
setMethod("Arith",signature(e1="real_symmetric_matrix" ,e2="numeric"              ),    rsm_arith_numeric)
setMethod("Arith",signature(e1="numeric"               ,e2="real_symmetric_matrix"),numeric_arith_rsm    )

setMethod("[",signature(x="real_symmetric_matrix",i="index",j="missing",drop="logical"),
          function(x,i,j,drop){
              out <- as.matrix(x)[,i,drop=FALSE]
              if(drop){
                  if(ncol(out)==1){
                      return(vec_to_rsm1(c(out)))
                  } else {
                      stop("for >1 element, use as.list()")
                  } 
              } else {
                  return(as.jordan(out,x))
              }
          } )
              

setReplaceMethod("[",signature(x="real_symmetric_matrix",i="index",j="missing",value="numeric"),function(x,i,j,value){stop("not defined")})
setReplaceMethod("[",signature(x="real_symmetric_matrix",i="index",j="missing",value="real_symmetric_matrix"),
                 function(x,i,j,value){
                   out <- as.matrix(x)
                   out[,i] <- as.matrix(value)  # the meat
                   return(as.jordan(out,x))
                 } )

setReplaceMethod("[",signature(x="real_symmetric_matrix",i="index",j="ANY",value="ANY"),function(x,i,j,value){stop("second argument redundant")})
