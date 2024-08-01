#include "Motif.hpp"

MotifSobol::MotifSobol(bool transformed): MotifPure(), _transformed(transformed){};

KMA::matrix MotifSobol::compute_v_new(const KMA::Mfield& Y_inters_k,
                                      const KMA::umatrix & Y_inters_supp,
                                      const arma::urowvec & v_dom,
                                      arma::uword v_len,
                                      const KMA::vector & p_k,
                                      arma::uword d,
                                      arma::uword m) const
{
  KMA::matrix v_new(v_len,d,arma::fill::zeros);
  if (Y_inters_k.n_rows == 1){
    v_new.rows(arma::find(v_dom==1)) = Y_inters_k(0);
    v_new.rows(arma::find(v_dom==0)).fill(arma::datum::nan);
    return v_new;
  } 
  KMA::vector coeff_k = arma::pow(p_k.elem(arma::find(p_k > 0)), m) / arma::sum(Y_inters_supp, 1);
  KMA::vector coeff_x = arma::sum(arma::conv_to<KMA::matrix>::from(Y_inters_supp) % arma::repmat(coeff_k,1,Y_inters_supp.n_cols), 0).t();
  coeff_x.elem(arma::find(arma::sum(Y_inters_supp, 0) == 0)).fill(arma::datum::nan);
  for (arma::uword i = 0; i < Y_inters_k.n_rows; ++i){
    v_new.rows(arma::find(v_dom==1)) += (Y_inters_k(i)*coeff_k(i)) / arma::repmat(coeff_x,1,Y_inters_k(i).n_cols);
  }
  
  v_new.rows(arma::find(v_dom==0)).fill(arma::datum::nan);
  return v_new;
}





