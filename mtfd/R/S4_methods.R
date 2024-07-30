#' @title showDetails
#' @description Shows details of an object of class MyS4Class.
#' @param object An object of class MyS4Class.
#' @export
setGeneric("generateCurves", function(object,sd_noise)
  standardGeneric("generateCurves")
)

#' @rdname showDetails
#' @export
setMethod("generateCurves", "motifSimulation", function(object,sd_noise) {
  breaks=lapply(rep(object@len,length.out=object@N),seq,from=0,by=object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
  basis=lapply(breaks,function(breaks_i) create.bspline.basis(norder=object@norder,breaks=breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
  len_motifs <- unlist(lapply(object@mot_details,function(x){x$len}))
  # loop for each curve
  fd_curves <- NULL
  if(!all(sapply(object@motifs_in_curves, is.null))) {
    fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,len_motifs){
      coeff=rep(NA,len_i) # coefficients of the curve = degree of freedom are initialized null
      if(object@distribution=='unif'){
        if(is.null(motifs_in_curves_i)){
          coeff=runif(len_i,min=object@coeff_min,max=object@coeff_max) # If we don't have motifs we sample len_i coefficients uniformely
        }else{
          pos_coeff_motifs=unlist(mapply(
            function(a,b) seq(a)+b,
            rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
            motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
          coeff[pos_coeff_motifs]= unlist(lapply(motifs_in_curves_i$motif_ids, function(id) {
            object@mot_details[[id]]$weights})) +
            rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
          # For each chosen coefficient add a uniform number and a gaussian noise
          coeff[-pos_coeff_motifs]=runif(len_i-length(pos_coeff_motifs),min=object@coeff_min,max=object@coeff_max) # All other curve coefficients are randomly generated
        }
        # Same as before but sampling from a beta
      }else if(object@distribution=='beta'){
        if(is.null(motifs_in_curves_i)){
          coeff=object@coeff_min+rbeta(len_i,0.45,0.45)*(object@coeff_max-object@coeff_min)
        }else{
          pos_coeff_motifs=unlist(mapply(
            function(a,b) seq(a)+b,
            rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
            motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
          coeff[pos_coeff_motifs]=unlist(lapply(motifs_in_curves_i$motif_ids, function(id) {
            object@mot_details[[id]]$weights})) +
            rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_ids]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
          # For each chosen coefficient add a uniform number and a gaussian noise
          coeff[-pos_coeff_motifs]=runif(len_i-length(pos_coeff_motifs),min=object@coeff_min,max=object@coeff_max) # All other curve coefficients are randomly generated
        }
      }else{
        stop('Wrong \'distrib\'')
      }
      curve=fd(coef=coeff,basisobj=basis_i) # Fitting curves using such coefficients and basis
      return(curve)
    },object@motifs_in_curves,object@norder-1+rep(object@len/object@dist_knots,length.out=object@N),basis,MoreArgs = list(len_motifs),SIMPLIFY=FALSE)
  }
  return(fd_curves=fd_curves)
})