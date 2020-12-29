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
         stop(paste("comparision operator \"", .Generic, "\" not defined for jordans"))
         )
}

`is.zero` <- function(e1,e2=0){
  stopifnot(is.numeric(e2))
  stopifnot(length(e2)==1)
  stopifnot(round(e2)==e2)
  stopifnot(e2==0)
  apply(as.matrix(e1),2,function(x){all(x==0)})
}

`jordan_compare_numeric` <- function(e1,e2){
   out <- is.zero(e1,e2)  # the meat

   switch(.Generic,
          "==" =  out,
          "!=" = !out,
          stop(paste("comparision operator \"", .Generic, "\" not defined for jordans"))
          )
}

`numeric_compare_jordan` <- function(e1,e2){
   out <- is.zero(e2,e1) # the meat; NB e1,e2 swapped WRT jordan_compare_numeric()

   switch(.Generic,
          "==" =  out,
          "!=" = !out,
          stop(paste("comparision operator \"", .Generic, "\" not defined for jordans"))
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

setMethod("show","jordan_matrix",
          function(object){
              object <- as.matrix(object)
              if(is.null(colnames(object))){
                  colnames(object) <- paste("[",seq_len(ncol(object)),"]",sep="")
                  }
              print(as.matrix(object))
              return(object)
          } )

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
