`quadraticform` <- function(M){ # modelled on lorentz::sol()
  if(missing(M)){  # return quadratic form
    jj <- getOption("quadraticform")
    if(is.null(jj)){
      cat("identity matrix\n")
      return(invisible(NULL))
    } else {
      return(jj)
    }
  } else { # set quadratic form; notionally a symmetric matrix
    stopifnot(is.matrix(M))
    stopifnot(nrow(M) == ncol(M))
    options("quadraticform" = M)
    return(M)
  }
}

`spin` <- function(a,V){
  stopifnot(is.numeric(a))
  stopifnot(is.matrix(V))
  out <- rbind(a,V)
  rownames(out) <- NULL
  class(out) <- "spin"   # this is the only place that class spin is assigned
  return(out)
}

`is.spin` <- function(x){inherits(x,"spin")}
`as.spin` <- function(x,d){
  if(is.spin(x)){
    return(x)
  } else if(is.matrix(x)){
    return(spin(x[1,,drop=TRUE],x[-1,,drop=FALSE]))
  } else if(is.numeric(x) & is.vector(x)){
    return(spin(x,matrix(0,d,length(x))))
  } else if(is.list(x)){
    return(spin(x[[1]],x[[2]]))
  } else {
    stop("not recognised")
  }
}

setGeneric("dim")
`dim.spin` <- function(x){ nrow(rn(x)) }

`length.spin` <- function(x){ncol(unclass(x))}

setGeneric("names")
`names.spin` <- function(x){colnames(unclass(x))}

setGeneric("names<-")
`names<-.spin` <- function(x,value){
  a <- r1(x)
  names(a) <- value
  return(spin(a,rn(x)))
}

`r1` <- function(x){unclass(x)[ 1,,drop=TRUE ]}
`rn` <- function(x){unclass(x)[-1,,drop=FALSE]}

`rspin` <- function(n=5,d=7,s=9){ spin(sample(s,n,replace=TRUE),matrix(sample(s,n*d,replace=TRUE),d,n))}

`print.spin` <- function(x){

  out <- rbind(r1(x),"--",rn(x))
  out[] <- format(c(out),width=max(nchar(c(out))),justify="right")
  jj <- c("1","",paste("[",seq_len(nrow(rn(x))),"]",sep=""))
  rownames(out) <- format(jj,width=max(nchar(jj)),justify="right")

  print(noquote(out))
  return(x)
}

`Ops.spin` <-
  function (e1, e2 = NULL) 
{
  f <- function(...){stop("odd---neither argument has class spin?")}
  unary <- nargs() == 1
  lclass <- nchar(.Method[1]) > 0
  rclass <- !unary && (nchar(.Method[2]) > 0)

  if(lclass & rclass){stopifnot(dim(e1) == dim(e2))}
  
  if(unary){
    if (.Generic == "+") {
      return(e1)
    } else if (.Generic == "-") {
      return(Sneg(e1))
    } else {
      stop("Unary operator '", .Generic, "' is not implemented for spin objects")
    }
  }
  if (!is.element(.Generic, c("+", "-", "*", "/", "^", "==", "!="))){
    stop("operator '", .Generic, "' is not implemented for spin objects")
  }

  if (.Generic == "*") {
    if (lclass && rclass) {
      return(SprodS(e1, e2))
    } else if (lclass) {
      return(SprodS(e1, as.spin(e2,dim(e1))))
    } else if (rclass) {
      return(SprodS(as.spin(e1,dim(e2)), e2))
    } else {
      f()
    }

  } else if (.Generic == "+") { 
    if (lclass && rclass) {
      return(SplusS(e1, e2))
    } else if (lclass) {
      return(SplusS(e1, as.spin(e2,dim(e1))))
    } else if (rclass) {
      return(SplusS(as.spin(e1,dim(e2)), e2))
    } else {
      f()
    }

  } else if (.Generic == "-") { 
    if (lclass && rclass) {
      return(SplusS(e1, Sneg(e2)))
    } else if (lclass) {
      return(SplusS(e1, as.spin(-e2,dim(e1))))
    } else if (rclass) {
      return(SplusS(as.spin(-e1,dim(e2)), e2))
    } else {
      f()
    }

  } else if (.Generic == "/") {
    if (lclass && rclass) {
      return(SdivS(e1,e2))    # error
    } else if (lclass) {
      return(SprodS(e1,as.spin(1/e2,dim(e1)))) # works 
    } else if (rclass) {
      return(SdivS(e1,e2))    # error
    } else {
      f()
    }
    
  } else if (.Generic == "^") {
    if (lclass && rclass) {
      return(SpowerS(e1,e2)) # error  
    } else if (lclass) {
      return(SpowerN(e1,e2)) # works
    } else if (rclass) {
      return(SpowerS(e1,e2)) # error
    } else {
      f()
    }

  } else if (.Generic == "==") {
    return(S.eq.S(e1,e2))
  } else if (.Generic == "!=") {
    return(!S.eq.S(e1,e2))
  } else {
    stop("should not reach here")
  }
}

`sync` <- function(e1,e2){
  jj <- rbind(seq_along(e1),seq_along(e2))
  e1 <- unclass(e1)
  e2 <- unclass(e2)
  e1 <- e1[,jj[1,],drop=FALSE]
  e2 <- e2[,jj[2,],drop=FALSE]
  list(
      s1 = e1[ 1,,drop=TRUE ],
      s2 = e2[ 1,,drop=TRUE ],
      v1 = e1[-1,,drop=FALSE],
      v2 = e2[-1,,drop=FALSE]
  )
}

`SprodS`  <- function(e1,e2){
  if(is.null(getOption("quadraticform"))){
    innerprod <- function(v1,v2){colSums(v1*v2)}
  } else {
    innerprod <- function(v1,v2){emulator::quad.3diag(quadraticform(),v1,v2)}
  }

  with(sync(e1,e2), spin(s1*s2+innerprod(v1,v2), sweep(v2,2,s1,"*")+sweep(v1,2,s2,"*")))
}

`Sneg`    <- function(e1){ spin(-r1(e1),-rn(e1))  }
`SplusS`  <- function(e1,e2){ with(sync(e1,e2), spin(s1+s2,v1+v2)) }
`S.eq.S`  <- function(e1,e2){  with(sync(e1,e2), (s1==s2) & apply(v2==v2,2,all)) }

`SdivS`   <- function(...){ stop("not a division algebra") }
`SpowerS` <- function(...){ stop("spin^spin not defined") }

`Spower_single_n` <- function(e1,n){
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

`SpowerN` <- function(e1,e2){
  n <- e2  # yes it's redundant but using e2 for n drives me nuts
  if(length(n)==1){
    return(Spower_single_n(e1,n))
  } else {
    jj <- rbind(seq_along(e1),seq_along(n))
    e1 <- unclass(e1)
    e1 <- e1[,jj[1,],drop=FALSE]
    n  <-  n[ jj[2,],drop=FALSE]
    for(i in seq_len(ncol(e1))){
      e1[,i] <- unclass(Spower_single_n(as.spin(e1[,i,drop=FALSE]),n[i]))
    }
    return(as.spin(e1))
  }
}

`[.spin` <- function(x, ...){ spin(r1(x)[...],rn(x)[,...,drop=FALSE]) }

`[<-.spin` <- function(x,index,value){
  x <- unclass(x)
  x[,index] <- unclass(as.spin(value,nrow(x)-1))
  as.spin(x)
}

