#ifndef __PerformanceIndex_ipp__
#define __PerformanceIndex_ipp__
#include "PerformanceIndex.hpp"

template<bool use1>
double PerformanceSobol::compute_Jk_helper(const KMA::Mfield& V,const KMA::ivector& s_k,
                                           const KMA::vector& p_k,const KMA::Mfield& Y,
                                           const KMA::vector& w,int m,double c_k, 
                                           const KMA::uvector & keep_k,
                                           const std::shared_ptr<Dissimilarity>& diss) const
{
  // domain of the centroid
  const arma::urowvec& v_dom = util::findDomain(V(0,0));
  
  // select the part of the domain of the centroid
  KMA::Mfield v_new = util::selectDomain<use1>(v_dom,V); 
  
  const unsigned int Y_size = Y.n_rows;
  
  arma::field<arma::mat> Y_inters_k(Y_size,1 + use1);

  shiftCurveHandle<use1>(Y_inters_k,Y,s_k,v_dom);
  
  if(std::isfinite(c_k) && (keep_k.size() > 0))
  { 
    arma::ivec supp_inters_length(arma::accu(keep_k));
    unsigned int c_k_ = static_cast<unsigned int>(c_k);

    unsigned int k = 0;
    for (arma::uword i = 0; i < Y_size; ++i){
      if (keep_k(i)){
        supp_inters_length(k++)= arma::accu(util::findDomain(Y_inters_k(i,0)));
      }
    }
    
    arma::uvec check_lengths = supp_inters_length < c_k_;
    
    if (arma::accu(check_lengths) > 0){ 
      return NA_REAL;
    }
  }

  arma::vec dist(Y_size);
  for (arma::uword i = 0; i < Y_size; ++i){
    dist(i) = diss->computeDissimilarity(Y_inters_k.row(i),v_new); 
  }
  arma::vec result = dist % pow(p_k,m);
  return arma::accu(result.elem(arma::find_finite(result)));
}


#endif // __PerformanceIndex_ipp__