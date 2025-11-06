#ifndef __TYPE_TRAITS__
#define __TYPE_TRAITS__
#include "RcppArmadillo.h"

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
