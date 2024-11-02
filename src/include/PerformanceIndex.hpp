#ifndef __PERFORMANCE_INDEX__
#define __PERFORMANCE_INDEX__
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include "Utilities.hpp"
#include "Dissimilarity.hpp"
#include <memory>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class PerformanceIndexAB
{
public:
  
  virtual double compute_Jk(const KMA::Mfield& v,
                            const KMA::ivector& s_k,
                            const KMA::vector& p_k,
                            const KMA::Mfield& Y,
                            const KMA::vector& w,
                            int m,
                            double c_k, // actually is an int
                            const KMA::uvector & keep_k, 
                            const std::shared_ptr<Dissimilarity>& diss) const = 0; 

  protected:
  
  template<bool use1>
  void shiftCurveHandle(KMA::Mfield& Y_inters_k,
                        const KMA::Mfield& Y,
                        const KMA::ivector& s_k,
                        const arma::urowvec& v_dom) const
  {
    std::size_t Y_size = Y.n_rows;
    std::size_t v_len = v_dom.size();
    KMA::uvector indeces_dom = arma::find(v_dom==0);
    KMA::uvector filtered_j;
    unsigned int y_len;
    KMA::ivector index;
    int s_k_i;
   
    for (unsigned int i = 0; i < Y_size; ++i){
      s_k_i = s_k[i];
      index = arma::regspace<arma::ivec>(1, v_len - std::max(0, 1-s_k_i))+std::max(1,s_k_i)-1;
      Y_inters_k(i,0).set_size(v_len,Y(0,0).n_cols);
      y_len = Y(i,0).n_rows;
      filtered_j = arma::find(index <= y_len);
      Y_inters_k(i,0).fill(arma::datum::nan);
      Y_inters_k(i,0).rows(std::max(0, 1-s_k_i), std::max(0, 1-s_k_i) + filtered_j.n_elem - 1) =  Y(i,0).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
      Y_inters_k(i,0).shed_rows(indeces_dom);
      if constexpr(use1){
        Y_inters_k(i,1).set_size(v_len,Y(0,0).n_cols);
        Y_inters_k(i,1).fill(arma::datum::nan);
        Y_inters_k(i,1).rows(std::max(0, 1-s_k_i), std::max(0, 1-s_k_i) + filtered_j.n_elem - 1) =  Y(i,1).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
        Y_inters_k(i,1).shed_rows(indeces_dom);
      }
    } 
  }

};

class PerformanceSobol: public PerformanceIndexAB
{
protected:
  
  template<bool use1>
  double compute_Jk_helper(const KMA::Mfield& v,const KMA::ivector& s_k,
                           const KMA::vector& p_k,const KMA::Mfield& Y,
                           const KMA::vector& w,int m,double c_k, 
                           const KMA::uvector& keep_k,
                           const std::shared_ptr<Dissimilarity>& diss) const;
};


class PerformanceL2 final: public PerformanceSobol
{
  
  double compute_Jk(const KMA::Mfield& v,const KMA::ivector& s_k,
                    const KMA::vector& p_k,const KMA::Mfield& Y,
                    const KMA::vector& w,int m,double c_k, 
                    const KMA::uvector& keep_k,
                    const std::shared_ptr<Dissimilarity>& diss) const override;
};


class PerformanceH1 final: public PerformanceSobol
{
  
  double compute_Jk(const KMA::Mfield& v,const KMA::ivector& s_k,
                    const KMA::vector& p_k,const KMA::Mfield& Y,
                    const KMA::vector& w,int m,double c_k,
                    const KMA::uvector& keep_k,
                    const std::shared_ptr<Dissimilarity>& diss) const override;
  
};

#include "PerformanceIndex.ipp"

#endif //__PERFORMANCE_INDEX__
