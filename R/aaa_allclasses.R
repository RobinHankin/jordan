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

setAs(from="jordan_matrix",to="matrix",def=function(from){from@x})
setGeneric("as.matrix")
setMethod("as.matrix",signature(x="jordan_matrix"),function(x){as(x,"matrix")})

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
