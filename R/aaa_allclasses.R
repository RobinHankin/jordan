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

`jordan_compare` <- function(e1,e2){
  stopifnot(is.jordan(e1) | is.jordan(e2))
  jj <- harmonize_oo(e1,e2)
  out <- apply(jj[[1]]==jj[[2]],2,all)

  switch(.Generic,
         "==" =  out,
         "!=" = !out,
         stop(paste("comparision operator \"", .Generic, "\" not defined for jordans"))
         )
}

setMethod("Compare",signature(e1 = "jordan" , e2="jordan" ), jordan_compare)
setMethod("Compare",signature(e1 = "jordan" , e2="numeric"), jordan_compare)
setMethod("Compare",signature(e1 = "numeric", e2="jordan" ), jordan_compare)

setMethod("[", signature("jordan",i="index",j="missing",drop="ANY"),function(x,i,j,drop){as.albert(as.matrix(x)[,i,drop=FALSE])})
setMethod("[", signature("jordan",i="index",j="ANY",drop="ANY"),function(x,i,j,drop){stop("second indexing argument not needed")})



`jordan_compare` <- function(e1,e2){
    stopifnot(is.jordan(e1) | is.jordan(e2))
    
    if(!is.jordan(e1)){e1 <- as.jordan(e1,e2)}
    if(!is.jordan(e2)){e2 <- as.jordan(e2,e1)}
    jj <- harmonize_oo(e1,e2)
    out <- apply(jj[[1]]==jj[[2]],2,all)
    
    switch(.Generic,
           "==" =  out,
           "!=" = !out,
           stop(paste("comparison operator \"", .Generic, "\" not defined for onions"))
           )
}

setMethod("Compare",signature(e1="jordan" ,e2="jordan" ), jordan_compare)
setMethod("Compare",signature(e1="jordan" ,e2="numeric"), jordan_compare)
setMethod("Compare",signature(e1="numeric",e2="jordan" ), jordan_compare)
