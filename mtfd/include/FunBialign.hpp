#ifndef __FUNBIALIGN_HPP__
#define __FUNBIALIGN_HPP__

#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <vector>
#include <ranges>
#include <algorithm>
#include "TypeTraits.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class FunBialign
{
public:
  FunBialign(const KMA::matrix& data,
             unsigned int portion_len);
  
  double fMSR_adj(void) const;
  arma::sp_mat createDistance(void) const;
  
  KMA::matrix _windowData;
  
private:
  
  KMA::vector Hscore(void) const;
  
  Rcpp::IntegerVector deleteV(const Rcpp::List& list_of_recommendations_ordered,
                     const Rcpp::List&  all_accolites,
                     const Rcpp::NumericVector& vec_of_scores_ordered);
  
  unsigned int _outrows;
  unsigned int _outcols;
  std::vector<unsigned int> _numerosity; // max number of splitted curves 
  mutable KMA::vector _scoreData;
  bool _isSingle;
};


# endif //__FUNBIALIGN_HPP__