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

`is.jordan` <- function(x){inherits(x,"jordan")}

setAs(from="jordan",to="matrix",def=function(from){from@x})
setAs(from="jordan_matrix",to="matrix",def=function(from){from@x})
setMethod("as.matrix","jordan",function(x){as(x,"matrix")})
setMethod("as.matrix","jordan_matrix",function(x){as(x,"matrix")})

setGeneric("length")
setMethod("length","jordan",function(x){ncol(as(x,"matrix"))})

setGeneric("names")
setMethod("names","jordan",function(x){colnames(as(x,"matrix"))})

setGeneric("names<-")
setReplaceMethod("names","jordan",
                 function(x,value){
                   jj <- as.matrix(x)
                   colnames(jj) <- value
                   return(albert(jj))
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
