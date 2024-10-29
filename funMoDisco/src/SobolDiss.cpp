#include "Dissimilarity.hpp"

SobolDiss::SobolDiss(const KMA::vector& w, bool transformed):Dissimilarity(),_w(w), _transformed(transformed){};

double SobolDiss::distance(const KMA::matrix& y,
                           const KMA::matrix& v) const
{
  KMA::matrix diff = arma::square(y - v); //(y-v)^2
  
  diff.replace(arma::datum::nan,0);
  
  const arma::rowvec & col_sum = arma::sum(diff,0); //colSums
  
  unsigned int n_rows = 0;
  for(arma::uword i = 0; i < y.n_rows; ++i)
    if(is_finite(y.row(i))){n_rows += 1;}
    
  arma::urowvec length_dom(y.n_cols,arma::fill::value(n_rows)); //length of the domain
    
  return sum((col_sum/length_dom)%_w.t())/y.n_cols;
}
