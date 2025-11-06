#ifndef __MOTIF_IPP__
#define __MOTIF_IPP__
#include "Motif.hpp"

template<bool use1>
std::variant<MotifPure::indexField,KMA::Mfield>
MotifSobol::compute_motif_helper(const arma::urowvec& v_dom,
                                 const KMA::ivector& s_k,
                                 const KMA::vector& p_k,
                                 const KMA::Mfield& Y,
                                 double m) const

{
  if(arma::accu(p_k) == 0) {
      Rcpp::stop("Motif with no members! Degenerate cluster!");
  }

  arma::uword v_len = v_dom.n_elem; //length of the current domain

  arma::uword d = Y(0,0).n_cols; //number of cols of the curves

  KMA::uvector p_k_pos = arma::find(p_k > 0); //we will consider only the curves assigned to motif k

  KMA::Mfield Y_inters_k(p_k_pos.n_elem,Y.n_cols); //list of shifted curves

  KMA::umatrix Y_inters_supp(p_k_pos.n_elem,v_len); // for each shifted curve contains an uvec with the domain

  KMA::uvector filtered_j;
  KMA::ivector index;
  KMA::uvector indeces_dom = arma::find(v_dom==0);
  unsigned int y_len;

  for (arma::uword i = 0; i < p_k_pos.n_elem; ++i){
    y_len = Y(p_k_pos(i),0).n_rows; //length of the curve
    index = std::max(1,s_k(p_k_pos(i))) - 1 + arma::regspace<KMA::ivector>(1,v_len - std::max(0,1 - s_k(p_k_pos(i))));
    filtered_j = arma::find(index <= y_len);
    Y_inters_k(i,0).set_size(v_len,d);
    Y_inters_k(i,0).fill(arma::datum::nan);
    Y_inters_k(i,0).rows(std::max(0,1 - s_k(p_k_pos(i))),std::max(0,1 - s_k(p_k_pos(i))) + filtered_j.n_elem - 1) =  Y(p_k_pos(i),0).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
    Y_inters_k(i,0).shed_rows(indeces_dom);
    Y_inters_supp.row(i) = util::findDomain(Y_inters_k(i,0));
    Y_inters_k(i,0).replace(arma::datum::nan,0);
    if constexpr(use1) {
      Y_inters_k(i,1).set_size(v_len,d);
      Y_inters_k(i,1).fill(arma::datum::nan);
      Y_inters_k(i,1).rows(std::max(0,1 - s_k(p_k_pos(i))),std::max(0,1 - s_k(p_k_pos(i))) + filtered_j.n_elem - 1) =  Y(p_k_pos(i),1).rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
      Y_inters_k(i,1).shed_rows(indeces_dom);
      Y_inters_k(i,1).replace(arma::datum::nan,0);
    }
    if (_transformed){
      Y_inters_k.row(i) = util::transform_curves<use1>(Y_inters_k.row(i));
    }
  }

  KMA::Mfield v_new(1,Y.n_cols);

  if constexpr (use1) {
    v_new(0,1) = compute_v_new(Y_inters_k.col(1),
                               Y_inters_supp,
                               v_dom,v_len,p_k,d,m);

  }

  v_new(0,0) = compute_v_new(Y_inters_k.col(0),
                             Y_inters_supp,
                             v_dom,v_len,p_k,d,m);


  KMA::uvector v_new_domain = arma::find(util::findDomain<KMA::matrix>(v_new(0,0)) == 1);
  arma::sword index_min = v_new_domain.min();
  arma::sword index_max = v_new_domain.max();

  v_new(0,0) = v_new(0,0).rows(index_min,index_max);

  if constexpr (use1){
      v_new(0,1) = v_new(0,1).rows(index_min,index_max);
  }
  if (index_min > 0) {
    Rcpp::Rcout<<"ATTENZIONE NO TESTATO MOTIF.IPP 78"<<std::endl;
    return std::make_pair(v_new,index_min);
  }

  return v_new;
}

template<bool use1>
void MotifSobol::elongation(KMA::Mfield& V_new,
                            std::vector<arma::urowvec> & V_dom,
                            KMA::imatrix& S_k,
                            const arma::vec & p_k,
                            const arma::ivec& len_elong_k,
                            const arma::uvec& keep_k,
                            double c,
                            const KMA::Mfield Y,
                            const unsigned int index,
                            const Parameters& param,
                            const std::shared_ptr<PerformanceIndexAB>& perf,
                            const std::shared_ptr<Dissimilarity>& diss) const
{
  if(len_elong_k.empty()) return;

  const KMA::Mfield& v_new_k = V_new.row(index);
  const arma::urowvec& v_dom_k = V_dom[index];
  const arma::ivec& s_k = S_k.col(index);

  // new vec with zero at the top
  arma::ivec len_elong_k_zero(len_elong_k.size() + 1, arma::fill::zeros);
  std::copy(len_elong_k.begin(), len_elong_k.end(), len_elong_k_zero.begin() + 1);

  // create a matrix whose column_i contains the vector s_k - len_elong_k_zero[i]
  arma::uword len_elong_k_zero_size = len_elong_k_zero.size();
  arma::imat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);

  for (arma::uword i=0; i < len_elong_k_zero_size;++i) {
    s_k_elong_left_right_temp.col(i) = s_k - len_elong_k_zero(i);
  }

  // create a sequence of integer from len_elong_k_zero.size() to 1
  arma::ivec reversedSequence = arma::regspace<arma::ivec>(len_elong_k_zero_size,-1,1);
  reversedSequence(0) -= 1;

  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence
  std::vector<arma::ivec> s_k_elong_left_right = util::repeat_elements(s_k_elong_left_right_temp, reversedSequence);

  std::vector<arma::ivec> len_elong_k_right_list(len_elong_k_zero_size);
  const int max_len_elong_k = len_elong_k.back();
  unsigned int v_dom_elong_size = 0;

  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i) {
    len_elong_k_right_list[i] = len_elong_k_zero.elem(find(len_elong_k_zero <=  max_len_elong_k - len_elong_k_zero(i)));
    v_dom_elong_size += len_elong_k_right_list[i].size();
  }

  //  v_dom_elong_left_right will be a vector of arma::uvec containing all the elongated version of v_dom_k
  std::vector<arma::urowvec> v_dom_elong_left_right(v_dom_elong_size);
  const int v_dom_k_len = v_dom_k.n_elem;
  unsigned int k = 0;
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i){
    const arma::ivec & leng_elong_right_vector = len_elong_k_right_list[i];
    const unsigned int leng_elong_left = len_elong_k_zero(i);
    for (unsigned int j = 0; j <  leng_elong_right_vector.size(); ++j) {
      arma::urowvec temp(leng_elong_left + v_dom_k_len + leng_elong_right_vector(j), arma::fill::ones);
      temp.cols(leng_elong_left, leng_elong_left + v_dom_k.n_elem - 1) = v_dom_k;
      v_dom_elong_left_right[k++] = temp;
    }
  }

  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  arma::uvec not_start_with_NA(v_elong_left_right_size,arma::fill::zeros);
  KMA::Mfield v_elong_left_right(v_elong_left_right_size,Y.n_cols);


  for (int i = 0; i < v_elong_left_right_size; i++) {
    const auto & pair_motif_shift = compute_motif_helper<use1>(v_dom_elong_left_right[i+1], s_k_elong_left_right[i], p_k, Y, param._m);
    if (std::holds_alternative<KMA::Mfield>(pair_motif_shift)){
      v_elong_left_right.row(i) = *(std::get_if<KMA::Mfield>(&pair_motif_shift));
      not_start_with_NA(i) = 1;
    }
  }

// #if HAS_RANGES
// C++20 version using ranges
//auto not_NA_index = std::views::iota(0, v_elong_left_right_size)
//                  | std::views::filter([&not_start_with_NA](int index_j){
//                      return not_start_with_NA(index_j);
//                  });
//#else
// C++17 fallback using simple loop: build a vector of valid indices
std::vector<int> not_NA_index_vec;
not_NA_index_vec.reserve(v_elong_left_right_size);
for (int index_j = 0; index_j < v_elong_left_right_size; ++index_j) {
    if (not_start_with_NA(index_j)) {
        not_NA_index_vec.push_back(index_j);
    }
}
//#endif

// Allocate filtered containers
KMA::Mfield filtered_v_elong(arma::accu(not_start_with_NA), Y.n_cols);
std::vector<arma::ivec> filtered_s_k(arma::accu(not_start_with_NA));

int i = 0;

#if HAS_RANGES
for (const auto &index : not_NA_index) {
    filtered_v_elong.row(i) = v_elong_left_right.row(index);
    filtered_s_k[i] = s_k_elong_left_right[index];
    ++i;
}
#else
for (int idx : not_NA_index_vec) {
    filtered_v_elong.row(i) = v_elong_left_right.row(idx);
    filtered_s_k[i] = s_k_elong_left_right[idx];
    ++i;
}
#endif


  const unsigned int s_k_elong_left_right_size = filtered_s_k.size();

  // compute performance index before elongation
  double Jk_before = perf->compute_Jk(v_new_k,s_k,p_k,Y,param._w,param._m,
                                      arma::datum::nan, KMA::uvector({}),
                                      diss);

  // compute performance indexes for all possible elongations
  arma::vec Jk_after(s_k_elong_left_right_size);
  double c_i;

  for (arma::uword i = 0; i < s_k_elong_left_right_size; i++) {

    c_i = std::max<double>(floor(util::findDomain(filtered_v_elong(i,0)).n_elem*(1 - param._max_gap)),c);

    Jk_after[i] = perf->compute_Jk(filtered_v_elong.row(i),filtered_s_k[i],
                                   p_k,Y,param._w,param._m,c_i,
                                   keep_k,diss);
  }

  // find the best elongation in terms of perf. index
  arma::vec diff_perc = (Jk_after-Jk_before)/Jk_before;
  arma::uword best_elong = arma::index_min(diff_perc);

  // check that the min really exists
  bool elongate = false;
  if (best_elong < diff_perc.size() and arma::conv_to<arma::Col<double>>::from(arma::find_nan(diff_perc)).size() != diff_perc.size())
    elongate = diff_perc(best_elong) < param._deltaJK_elong;

  // evaluate if elongate or not
  if(elongate) {
    V_new.row(index) =  filtered_v_elong.row(best_elong);
    V_dom[index] =  v_dom_elong_left_right[best_elong + 1];
    S_k.col(index) = filtered_s_k[best_elong]; // stampare questo
  } else {
    return;
  }
}

template<bool use1>
void MotifSobol::elongate_motifs_helper(KMA::Mfield& V_new,
                                        std::vector<arma::urowvec>& V_dom,
                                        KMA::imatrix& S_k,const KMA::matrix& P_k,
                                        const KMA::Mfield& Y,const KMA::matrix& D,
                                        const Parameters& param,
                                        const std::shared_ptr<PerformanceIndexAB>& perf,
                                        const std::shared_ptr<Dissimilarity>& diss,
                                        const Rcpp::Function & quantile_func) const
{
#ifdef _OPENMP
      // number of threads 
      const unsigned int n_threads = param._n_threads;
#endif

	std::size_t V_dom_size = V_dom.size();

	arma::uvec len_dom(V_dom_size);
	arma::uvec bool_gaps(V_dom_size);
  	for(unsigned int i = 0; i < V_dom_size; ++i){
    		len_dom(i) = V_dom[i].size();
    		bool_gaps(i) = (arma::accu(V_dom[i]) < len_dom(i));
  	}

  	arma::uvec with_gaps = arma::find(bool_gaps == 1);
  	unsigned int with_gaps_size = with_gaps.size();

  	// @TODO: check this part for domains with NaN
  	if (with_gaps_size != 0){

  		KMA::Mfield V_filled(with_gaps_size,Y.n_cols);

  		std::vector<KMA::uvector> V_dom_filled(with_gaps_size);

		  arma::vec Jk_before(with_gaps_size);

      arma::vec Jk_after(with_gaps_size);

     	// fill the domains of the motifs with gaps and recompute the motifs with the filled domains
     	// and compute the perf.indexes before and after the filling
		  #ifdef _OPENMP
		  #pragma omp parallel for num_threads(n_threads)
		  #endif
      for (unsigned int i = 0; i < with_gaps_size; ++i){

        V_dom_filled[i] = arma::uvec(V_dom[with_gaps(i)].n_elem,arma::fill::ones);
        const auto & variant_motif = compute_motif_helper<use1>(V_dom_filled[i],  // errore
                                                                		S_k.col(with_gaps(i)),
                                                                		P_k.col(with_gaps(i)),
                                                                		Y,param._m);

        V_filled.row(i) = *(std::get_if<KMA::Mfield>(&variant_motif)); // Possibile errore se ritorna Nullptr

        Jk_before(i) = perf->compute_Jk(V_new.row(with_gaps(i)),
                                        		S_k.col(with_gaps(i)),
                                        		P_k.col(with_gaps(i)),
                                        		Y,param._w,param._m,
                                        		arma::datum::nan,arma::uvec({}),
                                        		diss);

        Jk_after(i) = perf->compute_Jk(V_filled.row(i),
                                      S_k.col(with_gaps(i)),
                                      P_k.col(with_gaps(i)),
                                      Y,param._w,param._m,
                                      arma::datum::nan,arma::uvec({}),
                                      diss);
     }

		 // if filling the domain improves the perf. index over a certain threshold replace the domain and the motifs with the filled one
     const arma::uvec& fill = arma::find((Jk_after-Jk_before)/Jk_before < param._deltaJK_elong);
     for(unsigned int i=0; i < fill.size(); ++i){
       V_dom[with_gaps(i)] = V_dom_filled[i];
       V_new.row(with_gaps(i)) = V_filled.row(i);
       }
    }

    std::vector<arma::ivec> len_elong(V_dom_size);
    for (unsigned int i = 0; i < V_dom_size; ++i){
      const unsigned int len_max_elong_i = std::min<unsigned int>(std::floor(len_dom[i]*param._max_elong),
                                                                  param._c_max[i] - len_dom[i]);
      if (len_max_elong_i == 0){
        len_elong[i] = KMA::ivector{};
      }
      else{
        len_elong[i] =  (len_max_elong_i <= param._trials_elong) ?
        arma::regspace<arma::ivec>(1, len_max_elong_i):
        round(arma::linspace<arma::ivec>(1, len_max_elong_i, param._trials_elong));
      }
    }

    double quant = Rcpp::as<double>(quantile_func(D,0.25));
    // keep will be a matrix whose value i,j will be D(i,j) < quantile(0)
    arma::umat keep = D < quant;
    // col-wise sum of the matrix keep
    const arma::uvec& col_sum_keep = (sum(keep, 0)).t();
    // vector of bool = true iff col_sum_keep[i]==0
    const arma::uvec& col_sum_keep_zero = (col_sum_keep==0);
    // empty_k stores the indexes of the col that have col_sum_keep = 0
    const arma::uvec& empty_k = find(col_sum_keep_zero);

    for (auto k : empty_k){
      const unsigned int min_col_k_D = index_min(D.col(k));
      keep(min_col_k_D, k) = true;
    }

    #ifdef _OPENMP
    	#pragma omp parallel for num_threads(n_threads)
    #endif
    for (unsigned int i = 0; i < V_dom_size; ++i){
      elongation<use1>(V_new,V_dom,S_k,P_k.col(i),
                       len_elong[i],keep.col(i),param._c[i],
                       Y,i,param,perf,diss);
    }
}


#endif
