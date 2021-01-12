## quaternionic Hermitian matrices ("qhm"); setClass("quaternion_herm_matrix") in 'aaa_allclasses.R'

`quaternion_herm_matrix` <- function(M){new("quaternion_herm_matrix",x=cbind(M))}  # this is the only place new("real_symmetric_matrix",...) is called
`is.quaternion_herm_matrix` <- function(x){inherits(x,"quaternion_herm_matrix")}

`r_to_n_qhm` <- function(r){sqrt(1+8*r)/4}
`n_to_r_qhm` <- function(n){n*(2*n-1)}

`is_ok_qhm` <- function(r){ # 'r' = number of rows in [rowwise] matrix
    jj <- sqrt(1+8*r)
    if(jj == round(jj)){
        return((1+jj)/4) # size of nxn quaternion hermitian matrix
    } else {
        stop("not correct")
    }
}

`valid_qhm` <- function(object){
    x <- object@x
    if(!is.numeric(x)){
        return("not numeric")
    } else if(!is.matrix(x)){
        return("not a matrix")
    } else if(is_ok_qhm(nrow(x)) < 0){
        return("must have appropriate size")
    } else {
        return(TRUE)
    }
}

setValidity("quaternion_herm_matrix", valid_qhm)

`as.quaternion_herm_matrix` <- function(x,d,single=FALSE){  # single modelled on as.onion()
    if(is.quaternion_herm_matrix(x)){
        return(x)
    } else if(is.matrix(x)){
        return(quaternion_herm_matrix(x))
    } else if(is.vector(x)){
        if(single){
            return(quaternion_herm_matrix(x))
        } else {
            return(numeric_to_quaternion_herm_matrix(x,d)) # problem! we do not know how big it is
        }
    } else {
        stop("not recognised")
    }
}

`numeric_to_quaternion_herm_matrix` <- function(x,d){stop("no unique coercion")}

`rqhm` <- function(n=3,d=5){quaternion_herm_matrix(matrix(round(rnorm(n*(2*d^2-d)),2),ncol=n))}
`qhm_id` <- function(n,d){as.quaternion_herm_matrix(kronecker(qhm1_to_vec(diag(nrow=d)),t(rep(1,n))))}

`vec_to_qhm1` <- function(x){
   r <- length(x)
   n <- (1+sqrt(1+8*r))/4
   stopifnot(n==round(n))
   out <- onionmat(as.quaternion(0),n,n)
   out[upper.tri(out,FALSE)] <- as.quaternion(matrix(x[-seq_len(n)],nrow=4))
   out <- out + ht(out)
   diag(out) <- x[seq_len(n)]
   return(out)  # quaternion hermitian matrix
}

`qhm1_to_vec` <- function(M){
    ind <- upper.tri(M,FALSE)
    c(Re(diag(M)),
      t(cbind(Re(M[ind]),
              i (M[ind]),
              j (M[ind]),
              k (M[ind])
              ) ) ) 
}

`vec_qhmprod_vec` <- function(x,y){
    x <- vec_to_qhm1(x)
    y <- vec_to_qhm1(y)
    qhm1_to_vec((cprod(x,y)+cprod(y,x))/2)
}

setMethod("as.1matrix","quaternion_herm_matrix",function(x,drop=TRUE){
    out <- lapply(seq_along(x), function(i){x[i,drop=TRUE]})
    if((length(x)==1) & drop){out <- out[[1]]}
    return(out)
} )

`qhm_inverse` <- function(x){stop("inverses for qhm objects not implemented")}

`qhm_prod_qhm`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- vec_qhmprod_vec(jj[[1]][,i],jj[[2]][,i])
    }
    return(as.jordan(out,e1))
}

`qhm_power_numeric` <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- qhm1_to_vec(mymatrixpower_onion(vec_to_qhm1(jj[[1]][,i]),jj[[2]][i])) # the meat
    }
    return(as.jordan(out,e1))
}




`qhm_arith_qhm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_jordan(e1, e2),
         "-" = jordan_plus_jordan(e1,jordan_negative(e2)),
         "*" = qhm_prod_qhm(e1, e2),
         "/" = qhm_prod_qhm(e1, qhm_inverse(e2)), # fails
         "^" = stop("qhm^qhm not defined"),
         stop(paste("binary operator \"", .Generic, "\" not defined for qhm"))
         )
}

`qhm_arith_numeric` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e1, e2),  
         "-" = jordan_plus_numeric(e1,-e2),  
         "*" = jordan_prod_numeric(e1, e2),
         "/" = jordan_prod_numeric(e1, 1/e2),
         "^" = qhm_power_numeric(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for qhm"))
         )
}

`numeric_arith_qhm` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e2, e1),  
         "-" = jordan_plus_numeric(-e2,e1),  
         "*" = jordan_prod_numeric(e2, e1),
         "/" = jordan_prod_numeric(qhm_inverse(e2),e1),
         "^" = jordan_power_jordan(e2, e1),
         stop(paste("binary operator \"", .Generic, "\" not defined for qhm"))
         )
}

setMethod("Arith",signature(e1 = "quaternion_herm_matrix", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = jordan_negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on qhm objects"))
                   )
          } )

setMethod("Arith",signature(e1="quaternion_herm_matrix",e2="quaternion_herm_matrix"),    qhm_arith_qhm    )
setMethod("Arith",signature(e1="quaternion_herm_matrix",e2="numeric"               ),    qhm_arith_numeric)
setMethod("Arith",signature(e1="numeric"               ,e2="quaternion_herm_matrix"),numeric_arith_qhm    )

setMethod("[",signature(x="quaternion_herm_matrix",i="index",j="missing",drop="logical"),
          function(x,i,j,drop){
              out <- as.matrix(x)[,i,drop=FALSE]
              if(drop){
                  if(ncol(out)==1){
                      return(vec_to_qhm1(c(out)))
                  } else {
                      stop("for >1 element, use as.list()")
                  } 
              } else {
                  return(as.jordan(out,x))
              }
          } )
              
setReplaceMethod("[",signature(x="quaternion_herm_matrix",i="index",j="missing",value="quaternion_herm_matrix"),  # matches rsm equivalent
                 function(x,i,j,value){
                   out <- as.matrix(x)
                   out[,i] <- as.matrix(value)  # the meat
                   return(as.jordan(out,x))
                 } )
