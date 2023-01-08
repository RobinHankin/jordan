`spin` <- function(a,V){
  stopifnot(is.numeric(a))
  stopifnot(is.matrix(V))
  new("spin",x=rbind(a,V))  # this is the only place new("spin",...) is called
}

`r1` <- function(x){x@x[1,,drop=TRUE]}
`rn` <- function(x){x@x[-1,,drop=FALSE]}

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

`is.spin` <- function(x){inherits(x,"spin")}
`as.spin` <- function(x,d){
  if(is.spin(x)){
    return(x)
  } else if(is.matrix(x)){
    return(spin(a=x[1,,drop=TRUE],V=x[-1,,drop=FALSE]))
  } else if(is.numeric(x) & is.vector(x)){
    return(spin(a=x,V=matrix(0,d,length(x))))
  } else if(is.list(x)){
    return(spin(a=x[[1]],V=x[[2]]))
  } else {
    stop("not recognised")
  }
}

setGeneric("dim")
setMethod("dim","spin",function(x){ nrow(rn(x)) })

# names() defined for jordan objects
`rspin` <- function(n=3,d=5){spin(round(rnorm(n),2),matrix(round(rnorm(n*d),2),d,n))}

`spin_id` <- function(n=3,d=5){as.spin(rbind(1,matrix(0,d,n)))}

setMethod("show", "spin", function(object){spin_show(object)})

`spin_show` <- function(x){
     if(length(x)==0){
     cat("Null vector of", description(x,plural=TRUE),"\n")
     return((x))
  }
  cat("Vector of",description(x,plural=TRUE), "with entries\n")
  x <- as(x,"matrix")
  rownames(x) <- paste("[",seq_len(nrow(x))-1,"]",sep="")
  if(is.null(colnames(x))){
      colnames(x) <- paste("[",seq_len(ncol(x)),"]",sep="")
  }

  o <- getOption("head_and_tail")
  if(is.null(o)){o <- c(5,3)}
  if(length(o)==1){o <- c(o,o)}

  jj <- capture.output(x)
  n <- nrow(x)
  if(length(jj) > n+1){
      print(x)
      return(x)
  }

  substr(jj[2],1,3) <- " r "
  if(sum(o) < (n-1)){
      jj <- c(
          jj[1:2],
          paste(rep("-",nchar(jj[1])),collapse=""),
          jj[3:(o[1]+2)],
          paste(rep(".",nchar(jj[1])),collapse=""),
          jj[(n-o[2]+2):(n+1)]
      )
  }
  for(i in jj){
    cat(paste(i,"\n"))
  }
  return(x)
}

`harmonize_spin_spin` <- function(e1,e2){ # e1,e2: spin objects
  jj <- rbind(seq_along(e1),seq_along(e2))
  list(
      s1 = r1(e1)[ jj[1,]           ],
      s2 = r1(e2)[ jj[2,]           ],
      v1 = rn(e1)[,jj[1,],drop=FALSE],
      v2 = rn(e2)[,jj[2,],drop=FALSE]
  )
}

`harmonize_spin_numeric` <- function(e1,e2){ # e1: spin, e2: numeric
  jj <- rbind(seq_along(e1),seq_along(e2))
  list(
      s1 = r1(e1)[ jj[1,]],
      s2 =    e2 [ jj[2,]],
      v1 = rn(e1)[,jj[1,],drop=FALSE]
  )
}

`spin_prod_spin`  <- function(e1,e2){
  if(is.null(getOption("quadraticform"))){
    innerprod <- function(v1,v2){colSums(v1*v2)}
  } else {
    innerprod <- function(v1,v2){emulator::quad.3diag(quadraticform(),v1,v2)}
  }
  with(harmonize_spin_spin(e1,e2),{
    return(spin(a=s1*s2 + innerprod(v1,v2), V=sweep(v2,2,s1,"*")+sweep(v1,2,s2,"*")))})
}

`spin_prod_numeric`  <- function(e1,e2){with(harmonize_spin_numeric(e1,e2),{return(spin(a=s1*s2,V=sweep(v1,2,s2,"*")))})}

`spin_plus_numeric` <- function(e1,e2){stop("not implemented")}

`spin_negative`  <- function(e1){spin(-r1(e1),-rn(e1))}

`spin_plus_spin`  <- function(e1,e2){with(harmonize_spin_spin(e1,e2),{return(spin(s1+s2,v1+v2))})}

`spin_equal_spin`  <- function(e1,e2){with(harmonize_spin_spin(e1,e2),{return(spin(s1+s2,v1+v2))})}

`spin_inverse`   <- function(...){ stop("not a division algebra") }
`spin_power_spin` <- function(...){ stop("spin^spin not defined") }

`spin_power_single_n` <- function(e1,n){  # n a length-one nonnegative integer
  stopifnot(n==round(n))
  stopifnot(n>=0)
  stopifnot(length(n)==1)
  if(n==0){
    return(spin(a=1+0*r1(e1),V=rn(e1)*0)) # 1
  } else if(n==1){
    return(e1)
  } else { 
    ## return(e1*Recall(e1,n-1))  #  this would be inefficient
    out <- e1
    for(i in seq_len(n-1)){  # NB: n>=2
      out <- out*e1
    }
    return(out)
  }
}

`spin_power_numeric` <- function(e1,e2){
  stop("not yet implemented (it makes sense but I have not got round to implementing it yet)")
  n <- e2  # yes it's redundant but using e2 for n drives me nuts
  if(length(n)==1){
    return(spin_power_single_n(e1,n))
  } else {
    jj <- harmonize_spin_numeric(e1,n)

    }
    return(as.spin(e1))
  }

setMethod("Arith",signature(e1 = "spin", e2="missing"),
          function(e1,e2){
            switch(.Generic,
                   "+" = e1,
                   "-" = spin_negative(e1),
                   stop(gettextf("unary operator %s not defined for spin objects", dQuote(.Generic)))
                   )
          } )

## binary operators:

setMethod("Arith",signature(e1 = "spin", e2="spin"),
          function(e1,e2){
            switch(.Generic,
         "+" = spin_plus_spin(e1,  e2),
         "-" = spin_plus_spin(e1, spin_negative(e2)),
         "*" = spin_prod_spin(e1,  e2),
         "/" = stop("1/spin not defined"),
         "^" = stop("x^spin not defined"),
         stop(gettextf("binary operator %s not defined for spin objects", dQuote(.Generic)))
         )})

setMethod("Arith",signature(e1 = "spin", e2="numeric"),
          function(e1,e2){
            switch(.Generic,
         "+" = spin_plus_numeric(e1,e2),
         "-" = spin_plus_numeric(e1,-e2),
         "*" = spin_prod_numeric(e1,e2),
         "/" = spin_prod_numeric(e1,1/e2),
         "^" = spin_power_numeric(e1,  e2),
         stop(gettextf("binary operator %s not defined for spin objects", dQuote(.Generic)))
         )})

setMethod("Arith",signature(e1 = "numeric", e2="spin"),
          function(e1,e2){
            switch(.Generic,
         "+" = spin_plus_numeric(e2,e1),
         "-" = spin_plus_numeric(spin_negative(e2),e1),
         "*" = spin_prod_numeric(e2,e1),
         "/" = stop("1/spin not defined"),
         "^" = stop("x^spin not defined"),
         stop(gettextf("binary operator %s not defined for spin objects", dQuote(.Generic)))
         )})


setMethod("[", signature("spin",i="index",j="missing"),function(x,i,j,drop){spin(a=r1(x)[i],V=rn(x)[,i,drop=FALSE])})
setMethod("[", signature("spin",i="missing",j="index"),function(x,i,j,drop){stop()})
setMethod("[", signature("spin",i="missing",j="missing"),function(x,i,j,drop){x})
  
setReplaceMethod("[",signature(x="spin",i="index",j="missing",value="spin"),
                 function(x,i,j,value){
                   outa <- r1(x)
                   outa[i] <- r1(value)
                   outV <- rn(x)
                   outV[,i] <- rn(value)
                   return(spin(a=outa,V=outV))
                 } )

setReplaceMethod("[",signature(x="spin",i="index",j="missing",value="numeric"),
                 function(x,i,j,value){
                   stopifnot(length(value)==1)
                   stopifnot(value==0)
                   outa <- r1(x)
                   outa[i] <- value
                   outV <- rn(x)
                   outV[,i] <- 0 # the meat
                   return(spin(a=outa,V=outV))
                 } )

setReplaceMethod("[",signature(x="spin",i="ANY"  ,j="missing",value="ANY"),function(x,i,j,value){stop()})
setReplaceMethod("[",signature(x="spin",i="index"  ,j="index"  ,value="ANY"),function(x,i,j,value){stop()})
setReplaceMethod("[",signature(x="spin",i="missing",j="ANY"    ,value="numeric"),function(x,i,j,value){stop()})
setReplaceMethod("[",signature(x="spin",i="missing",j="missing",value="spin"),function(x,i,j,value){

  })
setReplaceMethod("[",signature(x="spin",i="missing",j="missing",value="numeric"),function(x,i,j,value){

  })
