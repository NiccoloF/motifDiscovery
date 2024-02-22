#include "Dissimilarity.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

L2::L2(const KMA::vector& w):SobolDiss(w) {};


double L2::computeDissimilarity(const KMA::Mfield& Y_i,
                                const KMA::Mfield& V_i) const
{
    return this->distance(Y_i(0,0),V_i(0,0));
}

void L2::set_parameters(const Parameters & newParameters){
    _w = newParameters._w;
}

void L2::computeDissimilarityClean(KMA::matrix & D_clean,
                                   const KMA::imatrix & S_clean,
                                   const std::vector<arma::urowvec> & V_dom_new,
                                   const KMA::Mfield & V_clean,
                                   const KMA::Mfield & Y) const
{
  return computeDissimilarityClean_helper<false>(D_clean,S_clean,V_dom_new,V_clean,Y);
}


KMA::vector L2::find_diss(const KMA::Mfield Y,
                          const KMA::Mfield V,
                          const KMA::vector& w,
                          double alpha, unsigned int c_k) const
{
  return find_diss_helper<false>(Y,V,w,alpha,c_k);
}
