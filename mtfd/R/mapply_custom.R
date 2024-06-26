#' @title .mapply_custom
#'
#' @description if the number of cluster is not null apply clusterMap, otherwise apply the classical version
#'
#' @author Marzia Angela Cremona & Francesca Chiaromonte
.mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    parallel::clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}