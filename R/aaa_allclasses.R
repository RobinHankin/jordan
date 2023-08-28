setClass("jordan", representation = "VIRTUAL" )

setClass("spin",
         slots    = c(x="matrix"),
         contains = "jordan"
         )

setClass("real_symmetric_matrix",
         slots = c(x="matrix"),
         contains = "jordan"
         )

setClass("complex_herm_matrix",
         slots = c(x="matrix"),
         contains = "jordan"
         )

setClass("quaternion_herm_matrix",
         slots = c(x="matrix"),
         contains = "jordan"
         )

setClass("albert",
         slots    = c(x="matrix"),  # a matrix with 27 rows
         contains = "jordan"
         )

setClassUnion("jordan_matrix", # everything except spin
              c("real_symmetric_matrix", "complex_herm_matrix",
                "quaternion_herm_matrix", "albert"))

setClassUnion("jordan_special", # everything except albert
              c("spin","real_symmetric_matrix", "complex_herm_matrix",
                "quaternion_herm_matrix"))

# setClassUnion("index", members =  c("numeric", "logical", "character")) # taken from the Matrix package.

`is.jordan` <- function(x){is(x,"jordan")}
`as.jordan` <- function(x,class){
    if(missing(class) & is.jordan(x)){return(x)}
    if(is.jordan(class)){class <- as.character(class(class))}
    switch(class,
           real_symmetric_matrix = as.real_symmetric_matrix(x),
           complex_herm_matrix = as.complex_herm_matrix(x),
           quaternion_herm_matrix = as.quaternion_herm_matrix(x),
           albert = as.albert(x),
           spin = as.spin(x),
           stop("not recognised")
           )
}

`as.identity` <- function(x){
           if(is.real_symmetric_matrix(x)){
               return(as.jordan(kronecker(   rsm1_to_vec(diag(nrow=nrow(x[1,drop=TRUE]))),t(rep(1,length(x)))),x))
           } else if(is.complex_herm_matrix(x)){
               return(as.jordan(kronecker(   chm1_to_vec(diag(nrow=nrow(x[1,drop=TRUE]))),t(rep(1,length(x)))),x))
           } else if(is.quaternion_herm_matrix(x)){
               return(as.jordan(kronecker(   qhm1_to_vec(diag(nrow=nrow(x[1,drop=TRUE]))),t(rep(1,length(x)))),x))
           } else if(is.albert(x)){
               return(as.jordan(kronecker(albert1_to_vec(diag(nrow=nrow(x[1,drop=TRUE]))),t(rep(1,length(x)))),x))
           } else if(is.spin(x)){
               return(spin(a=1+0*r1(x),V=rn(x)*0))
           } else {
               stop("not recognised")
           }
}

setAs(from="jordan",to="matrix",def=function(from){from@x})
setGeneric("as.matrix")
setMethod("as.matrix",signature(x="jordan"),function(x){as(x,"matrix")})

setGeneric("length")
setMethod("length","jordan",function(x){ncol(as.matrix(x))})

setGeneric("names")
setMethod("names","jordan",function(x){colnames(as.matrix(x))})

setGeneric("names<-")
setReplaceMethod("names","jordan",
                 function(x,value){
                   jj <- as.matrix(x)
                   colnames(jj) <- value
                   return(as.jordan(jj,as.character(class(x))))
                 } )

`jordan_compare_jordan` <- function(e1,e2){
  stopifnot(is.jordan(e1) | is.jordan(e2))
  jj <- harmonize_oo(e1,e2)
  out <- apply(jj[[1]]==jj[[2]],2,all)

  switch(.Generic,
         "==" =  out,
         "!=" = !out,
         stop(gettextf("comparison operator %s not defined for albert objects", dQuote(.Generic)))
         )
}


setGeneric("is.zero",function(x){standardGeneric("is.zero")})
setMethod("is.zero",signature(x="jordan"),function(x){is_zero_jordan(x)})

`is_zero_jordan` <- function(e1,e2=0){
  stopifnot(is.numeric(e2))
  stopifnot(length(e2)==1)
  stopifnot(round(e2)==e2)
  stopifnot(e2==0)
  apply(as.matrix(e1),2,function(x){all(x==0)})
}

`jordan_compare_numeric` <- function(e1,e2){
   out <- is_zero_jordan(e1,e2)  # the meat

   switch(.Generic,
          "==" =  out,
          "!=" = !out,
          stop(gettextf("comparison operator %s not defined for jordan objects", dQuote(.Generic)))
          )
}

`numeric_compare_jordan` <- function(e1,e2){
   out <- is_zero_jordan(e2,e1) # the meat; NB e1,e2 swapped WRT jordan_compare_numeric()

   switch(.Generic,
          "==" =  out,
          "!=" = !out,
          stop(gettextf("comparison operator %s not defined for jordan objects", dQuote(.Generic)))
          )
}

setMethod("Compare",signature(e1 = "jordan" , e2="jordan" ), jordan_compare_jordan)
setMethod("Compare",signature(e1 = "jordan" , e2="numeric"), jordan_compare_numeric)
setMethod("Compare",signature(e1 = "numeric", e2="jordan" ), numeric_compare_jordan)

setMethod("[", signature("jordan",i="index",j="missing",drop="ANY"),function(x,i,j,drop){as.jordan(as.matrix(x)[,i,drop=FALSE],x)})
setMethod("[", signature("jordan",i="index",j="ANY",drop="ANY"),function(x,i,j,drop){stop("second indexing argument not needed")})

## unary operators:
`jordan_negative` <- function(z){as.jordan(-as.matrix(z),z)}


## binary operators:
`jordan_plus_jordan`  <- function(e1,e2){
    jj <- harmonize_oo(e1,e2)
    as.jordan(jj[[1]] + jj[[2]],e1)
}

`jordan_plus_numeric`  <- function(e1,e2){
    jj <- harmonize_on(e1,e2)
    as.jordan(sweep(jj[[1]],2,jj[[2]],"+"),e1)
}

`jordan_prod_numeric` <- function(e1,e2){
    jj <- harmonize_on(e1,e2)
    as.jordan(sweep(jj[[1]],2,jj[[2]],"*"),e1)
}

setGeneric("length")
setMethod("length","jordan",function(x){ncol(as.matrix(x))})

setGeneric("sum")
setMethod("sum","jordan",function(x,na.rm=FALSE){as.jordan(cbind(rowSums(as.matrix(x))),x)})
   
setGeneric("as.1matrix",function(x,...){x})

setReplaceMethod("[",signature(x="jordan_matrix",i="index",j="missing",value="numeric"),
                 function(x,i,j,value){
                     stopifnot(length(value)==1)
                     stopifnot(value==0)
                     out <- as.matrix(x)
                     out[,i] <-  0  # the meat
                     return(as.jordan(out,x))
                 } )

`jordan_power_jordan` <- function(e1,e2){stop("x^jordan not defined")}

`mymatrixpower` <- function(M,n){
    stopifnot(length(n)==1)
    if(n==0){
        return(diag(nrow=nrow(M)))
    } else if(n==1){
        return(M)
    } else if(n<0){
        return(Recall(solve(M),-n))
    } else {
        jj <- eigen(M)
        D  <- jj$values
        stopifnot(is.numeric(D)) # verifies M is Hermitian
        O <- jj$vectors
        return(emulator::quad.tform(diag(D)^n,O))
    }
}

`mymatrixpower_onion` <- function(M,n){
    stopifnot(length(n)==1)
    if(n==0){
        return(M[1,1]*0 + diag(nrow=nrow(M)))
    } else if(n==1){
        return(M)
    } else if(n<0){
        stop("onion matrix inverses not implemented")
    } else {  # n>1
        out <- M
        for(i in seq_len(n-1)){ out <- out %*% M }  # the meat
        return(out)
    }
}

`description` <- function(x,plural=FALSE){
  if(is.real_symmetric_matrix(x)){
    out <- ifelse(plural,"real symmetric matrices","real symmetric matrix")
  } else if (is.complex_herm_matrix(x)){
    out <- ifelse(plural,"complex Hermitian matrices","complex Hermitian matrix")
  } else if (is.quaternion_herm_matrix(x)){
    out <- ifelse(plural,"quaternionic Hermitian matrices","quaternionic Hermitian matrix")
  } else if (is.albert(x)){
    out <- ifelse(plural,"Albert matrices","Albert matrix")
  } else if (is.spin(x)){
    out <- ifelse(plural,"spin objects","spin object")
  }  else {
    stop("not recognised")
  }
  return(out)
}

setMethod("show", "jordan_matrix", function(object){jordan_matrix_show(object)})

`jordan_matrix_show` <- function(x){
 if(length(x)==0){
     cat("Null vector of",description(x,plural=TRUE),"\n")
     return((x))
  }
  cat("Vector of",description(x,plural=TRUE), "with entries\n")
  x <- as(x,"matrix")
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
  } else {
      if(sum(o) < n){
          jj <- c(jj[seq_len(o[1]+1)],paste(rep(".",nchar(jj[1])),collapse=""),jj[(n-o[2]):(n+1)])
      }
      for(i in jj){
          cat(paste(i,"\n"))
      }
  }
 return(x)
}

conc_pair <- function(x,y){as.jordan(cbind(as.matrix(x),as.matrix(y)),x)}

setMethod("c","jordan",function(x,...){
  if (nargs() < 3){
    return(conc_pair(x, ...))
  } else {
    return(conc_pair(x, Recall(...)))
  }
} )

`matrix1_to_jordan` <- function(x){
  if(is.numeric(x)){
    return(as.real_symmetric_matrix(rsm1_to_vec(x),single=TRUE))
  } else if(is.complex(x)){
    return(as.complex_herm_matrix(chm1_to_vec(x),single=TRUE))
  } else if(is.onionmat(x)){
    jj <- x[1,1]
    if(is.quaternion(jj)){
      return(as.quaternion_herm_matrix(qhm1_to_vec(x),single=TRUE))
    } else if(is.octonion(jj)){
      return(as.albert(albert1_to_vec(x),single=TRUE))
    } else {
      stop("this cannot happen")
    }
  } else {
    stop("unrecognised matrix type")
  }
}
 

setGeneric("as.list")

setMethod("nrow","jordan",function(x){nrow(as.matrix(x))})
