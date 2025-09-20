#ifndef __TYPE_TRAITS__
#define __TYPE_TRAITS__
#include "RcppArmadillo.h"

#if __cpp_lib_ranges >= 201911L
#include <ranges>
#define HAS_RANGES 1
#else
#define HAS_RANGES 0
#endif

namespace KMA
{
  using matrix = arma::mat;
  using imatrix = arma::imat;
  using umatrix = arma::umat;
  
  using vector = arma::vec;
  using ivector = arma::ivec;
  using uvector = arma::uvec;
  
  using Vfield = arma::field<arma::vec>;
  using Mfield = arma::field<arma::mat>;
}


#endif // __TYPE_TRAITS__
