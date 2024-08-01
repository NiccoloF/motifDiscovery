#include "Dissimilarity.hpp"

template<bool use1>
KMA::vector SobolDiss::find_diss_helper(const KMA::Mfield Y,
                                        const KMA::Mfield V,
                                        const KMA::vector& w, 
                                        double alpha, unsigned int c_k) const
    {
      unsigned int d = Y(0,0).n_cols;
      arma::urowvec v_dom = util::findDomain(V(0,0));
      const KMA::Mfield& v_new = util::selectDomain<use1,arma::urowvec>(v_dom,V);
      int v_len = v_dom.size();
      int y_len = Y(0,0).n_rows;
      arma::ivec s_rep = arma::regspace<arma::ivec>(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
      const unsigned int s_rep_size = s_rep.size();
      KMA::Mfield y_rep(s_rep_size,Y.n_cols);
      
      KMA::uvector indeces_dom = arma::find(v_dom==0);
      KMA::ivector index;
      KMA::uvector filtered_j;
      KMA::uvector neg_index;
      
      for (unsigned int i = 0; i < s_rep_size; ++i) {
        index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
        filtered_j = arma::find((index > 0) && (index <= y_len));
        neg_index = arma::find(index <= 0); 
        y_rep(i,0).set_size(v_len, d);
        y_rep(i,0).fill(arma::datum::nan); 
        y_rep(i,0).rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = Y(0,0).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
        y_rep(i,0).shed_rows(indeces_dom); 
        if constexpr(use1) {
          y_rep(i,1).set_size(v_len, d);
          y_rep(i,1).fill(arma::datum::nan);
          y_rep(i,1).rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = Y(0,1).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
          y_rep(i,1).shed_rows(indeces_dom);
        }
      }
      
      KMA::ivector length_inter(s_rep_size);
      KMA::uvector non_na_indices;
      for (unsigned int i = 0; i < s_rep_size; ++i){
        non_na_indices = arma::find_nan(y_rep(i,0).col(0));
        length_inter(i) = y_rep(i,0).col(0).n_elem - non_na_indices.n_elem;
      }
      arma::uvec valid = length_inter >= c_k;
      if (arma::accu(valid) == 0) {        
        valid.elem(arma::find(length_inter == arma::max(length_inter))).fill(1);
      }
      
      double min_d = std::numeric_limits<double>::max();
      int min_s = 0;
    
      for (unsigned int i = 0; i < s_rep_size; i++) {
        if (valid(i)) {
        const double dist = computeDissimilarity(y_rep.row(i),v_new);
          if (dist < min_d){
            min_d = dist;
            min_s = s_rep[i];
          } 
        }
      }
      
      return arma::vec({static_cast<double>(min_s), min_d}); 
  }

template<bool use1>
KMA::vector SobolDiss::find_diss_aligned_helper(const KMA::Mfield Y,
                                                const KMA::Mfield V,
                                                bool aligned) const
{
      unsigned int d = Y(0,0).n_cols;
      arma::urowvec v_dom = util::findDomain(V(0,0));
      const KMA::Mfield& v_new = util::selectDomain<use1,arma::urowvec>(v_dom,V);
      int v_len = v_dom.size();
      int y_len = Y(0,0).n_rows;
      arma::ivec s_rep;
      if (aligned)
      {
        s_rep = 1;
      } 
      else 
      {
        s_rep = arma::regspace<arma::ivec>(1, y_len - v_len + 1);
      }
      const arma::uword s_rep_size = s_rep.size();
      KMA::Mfield y_rep(s_rep_size,Y.n_cols);
      
      KMA::uvector indeces_dom = arma::find(v_dom==0);
      KMA::ivector index;
      KMA::uvector filtered_j;
      KMA::uvector neg_index;
      for (arma::uword i = 0; i < s_rep_size; ++i) {
        index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
        filtered_j = arma::find((index > 0) && (index <= y_len));
        neg_index = arma::find(index <= 0); 
        y_rep(i,0).set_size(v_len, d);
        y_rep(i,0).fill(arma::datum::nan); 
        y_rep(i,0).rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = Y(0,0).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
        y_rep(i,0).shed_rows(indeces_dom); 
        if constexpr(use1) {
          y_rep(i,1).set_size(v_len, d);
          y_rep(i,1).fill(arma::datum::nan);
          y_rep(i,1).rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = Y(0,1).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
          y_rep(i,1).shed_rows(indeces_dom);
        }
      }
  
      double min_d = std::numeric_limits<double>::max();
      int min_s = 0;
      double dist;
      for (unsigned int i = 0; i < s_rep_size; i++) {
        dist = computeDissimilarity(y_rep.row(i),v_new);
        if (dist < min_d){
          min_d = dist;
          min_s = s_rep[i];
        } 
      }
      return arma::vec({static_cast<double>(min_s), min_d}); 
  }
  
template<bool use1>
void SobolDiss::computeDissimilarityClean_helper(KMA::matrix & D_clean,
                                                 const KMA::imatrix & S_clean,
                                                 const std::vector<arma::urowvec> & V_dom_new,
				                                         const KMA::Mfield & V_clean,										
                                                 const KMA::Mfield & Y) const
{ 
 const arma::uword d = Y(0,0).n_cols;
 KMA::Mfield y(1,Y.n_cols);
 arma::uword _n_rows_V = V_dom_new.size();
 arma::uword _n_rows_Y = Y.n_rows;

 for(arma::uword k=0; k < _n_rows_V; ++k)
 {
  const auto& s_k = S_clean.col(k);  
  const auto& v_dom = V_dom_new[k];  
  const int v_len = v_dom.size(); 
  const KMA::uvector & indeces_dom = arma::find(v_dom==0);
  KMA::ivector index;
  int y_len;
  for (arma::uword i=0; i < _n_rows_Y; ++i)
  {
    const int s = s_k(i);
    index = std::max(1,s) - 1 + arma::regspace<arma::ivec>(1,v_len - std::max(0,1-s));
    y_len = Y(i,0).n_rows; 
    y(0,0).set_size(v_len,d);
    y(0,0).fill(arma::datum::nan);
    const arma::uvec & filtered_j = arma::find(index <= y_len);
    y(0,0).rows(std::max(0, 1-s),std::max(0, 1-s) +  filtered_j.n_elem - 1) =  Y(i,0).rows(index(*(filtered_j.cbegin())) - 1,index(*(filtered_j.cend() - 1)) - 1);
    y(0,0).shed_rows(indeces_dom);
    if constexpr(use1) {
      y(0,1).set_size(v_len,d);
      y(0,1).fill(arma::datum::nan);
      y(0,1).rows(std::max(0, 1-s),std::max(0, 1-s) +  filtered_j.n_elem - 1) =  Y(i,1).rows(index(*(filtered_j.cbegin())) - 1,index(*(filtered_j.cend() - 1)) - 1);
      y(0,1).shed_rows(indeces_dom);
    }
   D_clean(i,k) = computeDissimilarity(y,V_clean.row(k)); 
  }
 }
}
