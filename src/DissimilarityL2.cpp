#include "Dissimilarity.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

L2::L2(const KMA::vector& w, bool transformed):SobolDiss(w, transformed) {};

void L2::set_both(bool both){ _both = both;}

double L2::computeDissimilarity(const KMA::Mfield& Y_i,
                                const KMA::Mfield& V_i) const
{
    if(_transformed && !_both)
    {
      const KMA::Mfield & Y_i_transf = util::transform_curves<false>(Y_i);
      return this->distance(Y_i_transf(0,0),V_i(0,0));
    }
    if(_transformed)
    { 
      const KMA::Mfield & Y_i_transf = util::transform_curves<false>(Y_i);
      const KMA::Mfield & V_i_transf = util::transform_curves<false>(V_i);
      return this->distance(Y_i_transf(0,0),V_i_transf(0,0));
    }
    return this->distance(Y_i(0,0),V_i(0,0));
}

void L2::set_parameters(const Parameters & newParameters){
    _w = newParameters._w;
    _transformed = newParameters._transformed;
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

KMA::vector L2::find_diss_aligned(const KMA::Mfield Y,
                                  const KMA::Mfield V,
                                  bool aligned) const
{
  return find_diss_aligned_helper<false>(Y,V,aligned);
}
