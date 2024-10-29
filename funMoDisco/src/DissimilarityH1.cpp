#include "Dissimilarity.hpp"

H1::H1(const KMA::vector& w,double alpha, bool transformed): SobolDiss(w, transformed), _alpha(alpha){};

void H1::set_both(bool both){ _both = both;}

double H1::computeDissimilarity(const KMA::Mfield& Y_i,
                                const KMA::Mfield& V_i) const
{
    if(_transformed && !_both)
    {
      const KMA::Mfield & Y_i_transf = util::transform_curves<true>(Y_i);
      return (1-_alpha) * this -> distance(Y_i_transf(0,0),V_i(0,0)) +
              _alpha * this -> distance(Y_i_transf(0,1),V_i(0,1));
    }
    if(_transformed){
      const KMA::Mfield & Y_i_transf = util::transform_curves<true>(Y_i);
      const KMA::Mfield & V_i_transf = util::transform_curves<true>(V_i);
      return (1-_alpha) * this -> distance(Y_i_transf(0,0),V_i_transf(0,0)) +
              _alpha * this -> distance(Y_i_transf(0,1),V_i_transf(0,1));
    }
    return (1-_alpha) * this -> distance(Y_i(0,0),V_i(0,0)) +
            _alpha * this -> distance(Y_i(0,1),V_i(0,1));
}

void H1::set_parameters(const Parameters & newParameters){
    _w = newParameters._w;
    _alpha = newParameters._alpha;
    _transformed = newParameters._transformed;
}

void H1::computeDissimilarityClean(KMA::matrix & D_clean,
                                   const KMA::imatrix & S_clean,
                                   const std::vector<arma::urowvec> & V_dom_new,
                                   const KMA::Mfield & V_clean,
                                   const KMA::Mfield & Y) const
{
  return computeDissimilarityClean_helper<true>(D_clean,S_clean,V_dom_new,V_clean,Y);
}

KMA::vector H1::find_diss(const KMA::Mfield Y,
                          const KMA::Mfield V,
                          const KMA::vector& w,
                          double alpha, unsigned int c_k) const
{
  return find_diss_helper<true>(Y,V,w,alpha,c_k);
}

KMA::vector H1::find_diss_aligned(const KMA::Mfield Y,
                                  const KMA::Mfield V,
                                  bool aligned) const
{
  return find_diss_aligned_helper<true>(Y,V,aligned);
}

