# Albert algebras; setClass("albert") is in  aaa_allclasses.R

`albert` <- function(M){new("albert",x=cbind(M))}  # this is the only place new("albert",...) is called
`is.albert` <- function(x){inherits(x,"albert")}

`r_to_n_albert` <- function(r=27){3}
`n_to_r_albert` <- function(n=3){27}

`is_ok_albert` <- function(r){
    if(r==27){
        return(3)
    } else {
        stop("not correct")
    }
}

`valid_albert` <- function(object){
  x <- object@x
  if(!is.numeric(x)){
    return("not numeric")
  } else if(!is.matrix(x)){
    return("not a matrix")
  } else if(!is_ok_albert(nrow(x))){
    return("must have 27 rows")
  } else {
    return(TRUE)
  }
}

setValidity("albert", valid_albert)

`as.albert` <- function(x,single=FALSE){  # single modelled on as.onion()
  if(is.albert(x)){
    return(x)
  } else if(is.matrix(x)){
    return(albert(x))
  } else if(is.vector(x)){
    if(single){
      return(albert(x))
    } else {
      return(numeric_to_albert(x))
    }
  } else {
    stop("not recognised")
  }
}

`ralbert` <- function(n=3){albert(matrix(round(rnorm(n*27),2),ncol=n))}
`albert_id` <- function(n){as.albert(kronecker(albert1_to_vec(herm_onion_mat(rep(1,3),as.octonion(rep(0,3)))),t(rep(1,n))))}

setMethod("show", "albert", function(object){albert_show(object)})
`albert_show` <- function(x){
  cat("Vector of",description(x,plural=TRUE), "with entries\n")
  jj <- as(x,"matrix")
  rownames(jj) <-
    c("    d1","    d2","    d3",
      "Re(o1)"," i(o1)"," j(o1)"," k(o1)"," l(o1)","il(o1)","jl(o1)","kl(o1)",
      "Re(o2)"," i(o2)"," j(o2)"," k(o2)"," l(o2)","il(o2)","jl(o2)","kl(o2)",
      "Re(o3)"," i(o3)"," j(o3)"," k(o3)"," l(o3)","il(o3)","jl(o3)","kl(o3)"
      )
  print(jj)
  return(x)
}

`albert_inverse` <- function(e1){stop("inverses not implemented in Albert algebras")}
`numeric_to_albert` <- function(e1){stop("There is no unique way to coerce numeric to Albert")}
`albert_prod_albert`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- vec_albertprod_vec(jj[[1]][,i],jj[[2]][,i])
    }
    return(albert(out))
}

`albert_arith_albert` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_jordan(e1, e2),
         "-" = jordan_plus_jordan(e1,jordan_negative(e2)),
         "*" = albert_prod_albert(e1, e2),
         "/" = albert_prod_albert(e1, albert_inverse(e2)), # fails
         "^" = stop("albert^albert not defined"),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

`albert_arith_numeric` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e1, e2),  
         "-" = jordan_plus_numeric(e1,-e2),  
         "*" = jordan_prod_numeric(e1, e2),
         "/" = jordan_prod_numeric(e1, 1/e2),
         "^" = albert_power_numeric(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

`numeric_arith_albert` <- function(e1,e2){
  switch(.Generic,
         "+" = jordan_plus_numeric(e2, e1),  
         "-" = jordan_plus_numeric(-e2,e1),  
         "*" = jordan_prod_numeric(e2, e1),
         "/" = jordan_prod_numeric(e2, 1/e1),
         "^" = albert_power_albert(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

setMethod("Arith",signature(e1 = "albert", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = jordan_negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on alberts"))
                   )
          } )

setMethod("Arith",signature(e1="albert" ,e2="albert" ), albert_arith_albert )
setMethod("Arith",signature(e1="albert" ,e2="numeric"), albert_arith_numeric)
setMethod("Arith",signature(e1="numeric",e2="albert" ),numeric_arith_albert )

`albert_power_albert` <- function(...){ stop("albert^albert not defined") }

`albert_power_numeric` <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- albert1_to_vec(mymatrixpower_onion(vec_to_albert1(jj[[1]][,i]),jj[[2]][i])) # the meat
    }
    return(as.jordan(out,e1))
}

`albert_power_single_n` <- function(e1,n){
    stopifnot(is.albert(e1))
    stopifnot(n==round(n))
    stopifnot(n>=0)
    stopifnot(length(n)==1)
    if(n==0){
        return(1+e1*0)
    } else if(n==1){
        return(e1)
    } else { 
        ## return(e1*Recall(e1,e2-1))  # inefficient
        out <- e1
        for(i in seq_len(n-1)){  # NB: n>=2
            out <- out*e1  # slightly inefficient as it does to matrix multiplications
        }
        return(out)
    }
}

`vec_to_albert1` <- function(x){
    stopifnot(length(x)==27)
    herm_onion_mat(x[1:3], as.octonion(matrix(x[-(1:3)],8,3)))
}

setMethod("as.1matrix","albert",function(x,drop=TRUE){
    out <- apply(as.matrix(x),2,vec_to_albert1,simplify=FALSE)
    if((length(x)==1) & drop){out <- out[[1]]}
    return(out)
} )

`vec_albertprod_vec` <- function(x,y){
  xmat <- vec_to_albert1(x)
  ymat <- vec_to_albert1(y)
  albert1_to_vec(cprod(xmat,ymat) + cprod(ymat,xmat))/2 ## xmat %*% ymat + ymat %*% xmat)/2
}

`albert1_to_vec` <- function(H){
    H <- H + onion::O0
  c(
      Re(c(H[1,1],H[2,2],H[3,3])),
      as.matrix(H[upper.tri(getM(H))])  # onion::getM()
  )
}

setMethod("as.list","albert", function(x){apply(as.matrix(x),2,vec_to_albert1)})

setMethod("[",signature(x="albert",i="index",j="missing",drop="logical"),
          function(x,i,j,drop){
              out <- as.matrix(x)[,i,drop=FALSE]
              if(drop){
                  if(ncol(out)==1){
                      return(vec_to_albert1(c(out)))
                  } else {
                      stop("for >1 element, use as.list()")
                  } 
              } else {
                  return(as.albert(out))
              }
          } )
              

setReplaceMethod("[",signature(x="albert",i="index",j="missing",value="albert"),
                 function(x,i,j,value){
                   out <- as.matrix(x)
                   out[,i] <- as.matrix(value)  # the meat
                   return(as.albert(out))
                 } )

