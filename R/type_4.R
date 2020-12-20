 ## Albert algebras; setClass("albert") is in  aaa_allclasses.R

`albert` <- function(M){new("albert",x=cbind(M))}  # this is the only place new("albert",...) is called
`is.albert` <- function(x){inherits(x,"albert")}

`valid_albert` <- function(object){
  x <- object@x
  if(!is.numeric(x)){
    return("not numeric")
  } else if(!is.matrix(x)){
    return("not a matrix")
  } else if(nrow(x) != 27){
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
      return(scalars_to_albert(x))
    }
  } else {
    stop("not recognised")
  }
}

`ralbert` <- function(n=5){albert(matrix(round(rnorm(n*27),2),nrow=27))}

setMethod("show", "albert", function(object){albert_show(object)})
`albert_show` <- function(x){
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


## unary operators:
`albert_negative` <- function(z){albert(-as.matrix(z))}
`albert_inverse` <- function(z){stop("inverses not implemented")}

## binary operators:
`albert_plus_albert`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    return(albert(jj[[1]] + jj[[2]]))
}

`albert_prod_albert`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    out <- jj[[1]]*0
    for(i in seq_len(ncol(out))){
        out[,i] <- v27_albertprod_v27(jj[[1]][,i],jj[[2]][,i])
    }
    return(albert(out))
}

`albert_plus_numeric` <- function(e1,e2){stop("not defined")}
`albert_prod_numeric` <- function(e1,e2){
    jj <- harmonize_on(e1,e2)
    as.albert(sweep(jj[[1]],2,jj[[2]],"*"))
}

`albert_arith_albert` <- function(e1,e2){
  switch(.Generic,
         "+" = albert_plus_albert(e1, e2),
         "-" = albert_plus_albert(e1,albert_negative(e2)),
         "*" = albert_prod_albert(e1, e2),
         "/" = albert_prod_albert(e1, albert_inverse(e2)), # fails
         "^" = stop("albert^albert not defined"),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

`albert_arith_numeric` <- function(e1,e2){
  switch(.Generic,
         "+" = albert_plus_numeric(e1, e2),  # fails
         "-" = albert_plus_numeric(e1,-e2),  # fails
         "*" = albert_prod_numeric(e1, e2),
         "/" = albert_prod_numeric(e1, 1/e2),
         "^" = albert_power_numeric(e1, e2),
         stop(paste("binary operator \"", .Generic, "\" not defined for alberts"))
         )
}

setMethod("Arith",signature(e1 = "albert", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = albert_negative(e1),
                   stop(paste("Unary operator", .Generic,
                              "not allowed on alberts"))
                   )
          } )


setMethod("Arith",signature(e1 = "albert",e2="albert"),albert_arith_albert)
setMethod("Arith",signature(e1 = "albert",e2="numeric"),albert_arith_numeric)
setMethod("Arith",signature(e1 = "albert",e2="albert"),albert_arith_albert)



`AL.eq.AL`  <- function(e1,e2){apply(unclass(e1) == unclass(e2),2,all)}
`ALprodS` <- function(e1,e2){ albert(unclass(e1)*e2)}

`AldivAL`   <- function(...){ stop("albert algebra not a division algebra") }
`AlpowerAL` <- function(...){ stop("albert^albert not defined") }

`ALpower_single_n` <- function(e1,n){
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
      out <- out*e1
    }
    return(out)
  }
}

`ALpowerN` <- function(e1,e2){
  n <- e2  # yes it's redundant but using e2 for n drives me nuts
  if(length(n)==1){
    return(Spower_single_n(e1,n))
  } else {
    jj <- rbind(seq_along(e1),seq_along(n))
    e1 <- unclass(e1)
    e1 <- e1[,jj[1,],drop=FALSE]
    n  <-  n[ jj[2,],drop=FALSE]
    for(i in seq_len(ncol(e1))){
      e1[,i] <- unclass(Spower_single_n(as.albert(e1[,i,drop=FALSE]),n[i]))
    }
    return(as.albert(e1))
  }
}

`sum.albert` <- function(x,na.rm=FALSE){ as.albert(cbind(rowSums(unclass(x)))) }
   
`v27_to_albertmatrix` <- function(x){
    stopifnot(length(x)==27)
    herm_onion_mat(x[1:3], as.octonion(matrix(x[-(1:3)],8,3)))
}

`v27_albertprod_v27` <- function(x,y){
  xmat <- v27_to_albertmatrix(x)
  ymat <- v27_to_albertmatrix(y)
  albertmatrix_to_v27(cprod(xmat,ymat) + cprod(ymat,xmat))/2 ## xmat %*% ymat + ymat %*% xmat)/2
}

`albertmatrix_to_v27` <- function(H){
  c(
      Re(c(H[1,1],H[2,2],H[3,3])),
      as.matrix(H[upper.tri(getM(H))])  # onion::getM()
  )
}

`real_to_albertmatrix` <- function(x){ 
  stopifnot(length(x)==1)
  herm_onion_mat(rep(x,3),onions=as.octonion(rep(0,3)))
}

scalars_to_albert <- function(x){
  out <- matrix(0,27,length(x))
  for(i in seq_along(x)){
    out[,i] <- albertmatrix_to_v27(real_to_albertmatrix(x[i]))
  }
  return(as.albert(out))
}

setGeneric("as.list")
setMethod("as.list","albert", function(x){apply(as.matrix(x),2,v27_to_albertmatrix)})

setReplaceMethod("[",signature(x="albert",i="index",j="missing",value="numeric"),
                 function(x,i,j,value){
                     out <- as.matrix(x)
                     out[,i] <-  as.matrix(as.albert(value))  # the meat
                     return(as.albert(out))
                 } )

setReplaceMethod("[",signature(x="albert",i="index",j="missing",value="albert"),
                 function(x,i,j,value){
                   out <- as.matrix(x)
                   out[,i] <- as.matrix(value)  # the meat
                   return(as.albert(out))
                 } )
