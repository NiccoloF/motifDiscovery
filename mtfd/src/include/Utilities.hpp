#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include <ranges>
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

  // to be generalized fully using template programming
  inline arma::rowvec apply_rmNA(const arma::mat & y0, 
                                 std::function<double(const arma::vec &)> func) { 
    arma::rowvec result(y0.n_cols);
    for (arma::uword j = 0; j < y0.n_cols; ++j)
    {
      const arma::vec & col_j = y0.col(j);
      result(j) = func(col_j.elem(arma::find_finite(col_j)));
    }
    return result;

  }

  template<bool use1>
  KMA::Mfield transform_curves(const KMA::Mfield &y) {
    KMA::Mfield y_transformed(1,1);
    if constexpr(use1)
      y_transformed.set_size(1,2);
    arma::rowvec y_min = apply_rmNA(y(0,0), [](const arma::vec & v){ return arma::min(v); }); 
    arma::rowvec y_max = apply_rmNA(y(0,0), [](const arma::vec & v){ return arma::max(v); }); 
    arma::rowvec diff = y_max - y_min;
    y_transformed(0,0) = ((y(0,0)- arma::repmat(y_min,y(0,0).n_rows,1))/ arma::repmat(diff,y(0,0).n_rows,1));
    y_transformed(0,0).cols(arma::find(diff == 0)).eval().fill(0.5); 
    if constexpr(use1)
    {
      y_transformed(0,1) = y(0,1)/ arma::repmat(diff,y(0,1).n_rows,1); 
      y_transformed(0,1).cols(arma::find(diff == 0)).eval().fill(0.0); 
    }
    return y_transformed;
  }

  template<typename T>
  arma::Mat<T> combn2(const arma::Col<T> & y){  // Col<double> = vec, Col<uword> = uvec, Col<sword> = ivec
    int n = y.n_elem;
    KMA::uvector v(n,arma::fill::zeros);  
    v(0) = 1;
    v(1) = 1;
    arma::uword l = 0;
    arma::Mat<T> result(2,n*(n-1)/2, arma::fill::zeros); 
    arma::uword k;
    do {
      k = 0;
      auto filter_index = std::views::iota(0,n) 
        | std::views::filter([&v](int i){return v(i);});
      for(auto i: filter_index) 
        result(k++,l) = y(i); 
      l++;
    } while (std::prev_permutation(v.begin(), v.end())); 
    return result;
  }
  
  template<class T>
  T repLem(const T & v,
           const arma::ivec & times){
    T result(arma::accu(times));
    arma::uword k = 0;
    arma::uword times_size = times.size();
    for (arma::uword i = 0; i < times_size; ++i) 
      for (arma::ivec::value_type j = 0; j < times[i]; ++j)
        result(k++) = v(i);
    return result;
  }


}
  


#endif // __UTILITIES_HPP__
  
