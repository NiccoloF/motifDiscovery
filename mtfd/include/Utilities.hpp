#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

namespace util
{
template <typename T>
concept IsArmaVector = std::is_same_v<T, KMA::uvector> || 
                       std::is_same_v<T, arma::urowvec>;

  template<bool use1,IsArmaVector T>
  KMA::Mfield selectDomain(const T& v_dom,const KMA::Mfield& V)
  { 
    KMA::uvector dom = arma::find(v_dom==1);
    KMA::Mfield v(1,V.n_cols);
    v(0,0) = V(0,0).rows(dom);
    
    if constexpr(use1)
      v(0,1) = V(0,1).rows(dom);
    
    return v;
  
  }
  
 
  // returns a rowvector 
  template <typename MatType>
  arma::urowvec findDomain(const MatType& v) {
    arma::urowvec result(v.n_rows, arma::fill::zeros);
    for (arma::uword i = 0; i < v.n_rows; ++i) {
      const KMA::uvector& finite_row = arma::find_finite(v.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return result;
  }
 
 
  inline std::vector<arma::ivec> repeat_elements(const KMA::imatrix& A,const KMA::ivector & times) {
    arma::uword times_size = times.n_elem;
    std::vector<KMA::ivector> result((times_size*(times_size+1))/2 - 1);
    std::size_t i = 0;
    for(arma::uword j = 0;j < times_size;++j)
    {
      const KMA::imatrix& B = arma::repmat(A.col(j),1,times[j]);
      B.each_col([&result,&i](const KMA::ivector& v){result[i++] = v;});
    }
    return result;
  }
  
}
  


#endif // __UTILITIES_HPP__