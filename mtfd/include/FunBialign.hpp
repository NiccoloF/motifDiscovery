#ifndef __FUNBIALIGN_HPP__
#define __FUNBIALIGN_HPP__

#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <vector>
#include <ranges>
#include <algorithm>
#include <execution>
#include "TypeTraits.hpp"
#include "fastcluster.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class FunBialign
{
public:
  FunBialign(const KMA::matrix& data,
             unsigned int portion_len);
  
  double fMSR_adj(void) const;
  KMA::vector createDistance(void) const;
  
private:
  
  KMA::vector Hscore(void) const;
  
  unsigned int _outrows;
  unsigned int _outcols;
  KMA::matrix _windowData;
  std::vector<unsigned int> _numerosity; // max number of splitted curves 
  mutable KMA::vector _scoreData;
  bool _isSingle;
};









# endif //__FUNBIALIGN_HPP__