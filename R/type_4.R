 ## Albert algebras

`albert` <- function(M){
  M <- cbind(M)
  stopifnot(is.matrix(M))
  stopifnot(nrow(M) == 27)
  class(M) <- "albert"   # this is the only place class albert is assigned
  return(M)
}
`is.albert` <- function(x){inherits(x,"albert")}

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

`length.albert` <- function(x){ncol(unclass(x))}

setGeneric("names")
`names.albert` <- function(x){colnames(unclass(x))}

`names<-.albert` <- function(x,value){
  x <- unclass(x)
  colnames(x) <- value
  return(albert(x))
}

`ralbert` <- function(n=5){albert(matrix(round(rnorm(n*27),2),nrow=27))}

`print.albert` <- function(x){
  jj <- unclass(x)
  rownames(jj) <-
    c("    d1","    d2","    d3",
      "Re(o1)"," i(o1)"," j(o1)"," k(o1)"," l(o1)","il(o1)","jl(o1)","kl(o1)",
      "Re(o2)"," i(o2)"," j(o2)"," k(o2)"," l(o2)","il(o2)","jl(o2)","kl(o2)",
      "Re(o3)"," i(o3)"," j(o3)"," k(o3)"," l(o3)","il(o3)","jl(o3)","kl(o3)"
      )
  print(jj)
  return(x)
}

`Ops.albert` <-
  function (e1, e2 = NULL) 
{
  f <- function(...){stop("odd---neither argument has class albert?")}
  unary <- nargs() == 1
  lclass <- nchar(.Method[1]) > 0
  rclass <- !unary && (nchar(.Method[2]) > 0)

  if(unary){
    if (.Generic == "+") {
      return(e1)
    } else if (.Generic == "-") {
      return(ALneg(e1))
    } else {
      stop("Unary operator '", .Generic, "' is not implemented for albert objects")
    }
  }
  if (!is.element(.Generic, c("+", "-", "*", "/", "^", "==", "!="))){
    stop("operator '", .Generic, "' is not implemented for albert objects")
  }

  if (.Generic == "*") {
    if (lclass && rclass) {
      return(ALprodAL(e1, e2))
    } else if (lclass) {
      return(ALprodS(e1, e2))
    } else if (rclass) {
      return(ALprodS(e2, e1))
    } else {
      f()
    }

  } else if (.Generic == "+") { 
    if (lclass && rclass) {
      return(ALplusAL(e1, e2))
    } else if (lclass) {
      return(ALplusS(e1, e2))
    } else if (rclass) {
      return(ALplusS(e2, e1))
    } else {
      f()
    }

  } else if (.Generic == "-") { 
    if (lclass && rclass) {
      return(ALplusAL(e1, ALneg(e2)))
    } else if (lclass) {
      return(ALplusS(e1, -e2))
    } else if (rclass) {
      return(ALplusS(e2, -e1))
    } else {
      f()
    }

  } else if (.Generic == "/") {
    if (lclass && rclass) {
      return(ALdivAL(e1,e2))   # error
    } else if (lclass) {
      return(ALprodS(e1,1/e2)) # works
    } else if (rclass) {
      return(ALdivAL(e1,e2))   # error
    } else {
      f()
    }
    
  } else if (.Generic == "^") {
    if (lclass && rclass) {
      return(ALpowerAL(e1,e2)) # error  
    } else if (lclass) {
      return(ALpowerN(e1,e2)) # works
    } else if (rclass) {
      return(ALpowerAL(e1,e2)) # error
    } else {
      f()
    }

  } else if (.Generic == "==") {
    return(AL.eq.AL(e1,e2))
  } else if (.Generic == "!=") {
    return(!AL.eq.AL(e1,e2))
  } else {
    stop("should not reach here")
  }
}

syncAL <- function(e1,e2){
  jj <- rbind(seq_along(e1),seq_along(e2))
  e1 <- unclass(e1)
  e2 <- unclass(e2)
  list(
      x1=e1[,jj[1,],drop=FALSE],
      x2=e2[,jj[2,],drop=FALSE]
       )
}

`ALplusAL`  <- function(e1,e2){with(syncAL(e1,e2), albert(x1+x2))}

`ALprodAL`  <- function(e1,e2){
  with(syncAL(e1,e2), {
    out <- x1*0
    for(i in seq_len(ncol(out))){
      out[,i] <- v27_albertprod_v27(x1[,i],x2[,i])
    }
    return(albert(out))
  })
}

`ALneg`    <- function(e1){albert(-unclass(e1))}
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

`[.albert` <- function(x, ...){
  out <- unclass(x)
  out <- out[, ..., drop=FALSE]
  albert(out)
}

`[<-.albert` <- function(x,index,value){
  out <- unclass(x)
  if(is.vector(value)){
    value <- kronecker(t(value),c(1,rep(0,nrow(x)-1)))
    value[2,] <- value[1,]
    value[3,] <- value[1,]
  }
  out[,index] <- value
  albert(out)
}

`sum.albert` <- function(x,na.rm=FALSE){ as.albert(cbind(rowSums(unclass(x)))) }
   
`v27_to_albertmatrix` <- function(x){
  herm_onion_mat(x[1:3], as.octonion(matrix(x[-(1:3)],8,3)))
}

`v27_albertprod_v27` <- function(x,y){
  xmat <- v27_to_albertmatrix(x)
  ymat <- v27_to_albertmatrix(y)
  albertmatrix_to_v27(onionmatprod(xmat,ymat) + onionmatprod(ymat,xmat))/2 ## xmat %*% ymat + ymat %*% xmat)/2
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

scalars_to_albertmatrix <- function(x){
  out <- matrix(0,27,length(x))
  for(i in seq_along(x)){
    out[i,] <- real_to_albertmatrix(x[i])
  }
  return(as.albert(out))
}
         
