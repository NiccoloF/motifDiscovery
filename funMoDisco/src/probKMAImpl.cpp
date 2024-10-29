#include "ProbKMA.hpp"
#include "Factory.hpp"
#include <forward_list>
#include <limits>
#include <string_view>

class ProbKMA::_probKMAImp
{
public:
    enum class StructType {motifs, curves};

    _probKMAImp(const Rcpp::List& Y,
                const Rcpp::List& parameters,
                const KMA::matrix& P0,const KMA::imatrix& S0,
                const std::string_view diss, 
                const Rcpp::List& V_init): 
                _probKMAImp(Y,parameters,P0,S0,diss)
                {
                  init_motifs = true;
                  initial_motifs(V_init, diss); 
                }; 

    _probKMAImp(const Rcpp::List& Y,
                const Rcpp::List& parameters,
                const KMA::matrix& P0,const KMA::imatrix& S0,
                const std::string_view diss):
                _parameters(parameters),_P0(P0),_S0(S0)
               {
                    // Initialize c++ Data structure
                    Initialize(Y,diss);

                    // Create Dissimilarity factory
                    util::SharedFactory<Dissimilarity> dissfac;
                    dissfac.FactoryRegister<L2>("L2",_parameters._w,_parameters._transformed); 
                    dissfac.FactoryRegister<H1>("H1",_parameters._w,_parameters._alpha,_parameters._transformed); 

                    //Create Motif factory
                    util::SharedFactory<MotifPure> motfac;
                    motfac.FactoryRegister<MotifL2>("L2",_parameters._transformed);
                    motfac.FactoryRegister<MotifH1>("H1",_parameters._transformed);

                    //Create Performance factory
                    util::SharedFactory<PerformanceIndexAB> perfac;
                    perfac.FactoryRegister<PerformanceL2>("L2");
                    perfac.FactoryRegister<PerformanceH1>("H1");

                    //check whether diss is valid
                    _motfac = motfac.instantiate(diss); // copy elision
                    _dissfac = dissfac.instantiate(diss); // copy elision
                    _perfac = perfac.instantiate(diss); // copy elision
                    if(not(_motfac and _dissfac and _perfac))
                      Rcpp::stop("Invalid dissimilarity: Choose between L2, H1");
                }

    ~_probKMAImp() = default;

    void Initialize(const Rcpp::List& Y,std::string_view diss)
    {
      // Convert Rcpp Data Structure(List) into Armadillo data Structure(field)
      const Rcpp::List& Y0 = Y[0];
      const Rcpp::List& Y1 = Y[1];

      if (diss == "H1") {
        handleCaseH1<StructType::curves>(Y0, Y1);
      } else if (diss == "L2") {
        handleCaseL2<StructType::curves>(Y0, Y1);
      }
      reinit_motifs(_parameters._c,_Y.front().n_cols);
    }

      // Support function for the case "H1"
      template<StructType T>
      void handleCaseH1(const Rcpp::List& Y0, const Rcpp::List& Y1) {
        const arma::uword Y_size = Y0.size();
        if constexpr (T == StructType::curves)
        {
          _Y.set_size(Y_size, 2);
        }
        else
        {
          _V.set_size(Y_size, 2);
        }

        for (arma::uword i = 0; i < Y_size; i++) {
          if constexpr (T == StructType::curves)
          {
            _Y(i, 0) = Rcpp::as<KMA::matrix>(Y0[i]);
            _Y(i, 1) = Rcpp::as<KMA::matrix>(Y1[i]);
          }
          else
          {
            _V(i, 0) = Rcpp::as<KMA::matrix>(Y0[i]);
            _V(i, 1) = Rcpp::as<KMA::matrix>(Y1[i]);
          }
        }
      }

      // Support function for the case "L2"
      template<StructType T>
      void handleCaseL2(const Rcpp::List& Y0, const Rcpp::List& Y1) {
        if (!Rf_isNull(Y0[0])) {
          _isY1 = false;
          handleNonNullY<T>(Y0);
        } else {
          _isY0 = false;
          handleNonNullY<T>(Y1);
        }
      }

      // Support function for the case "L2"
      template<StructType T>
      void handleNonNullY(const Rcpp::List& Y) {
        const arma::uword Y_size = Y.size();
        if constexpr (T == StructType::curves)
        {
          _Y.set_size(Y_size, 1);
        }
        else
        {
          _V.set_size(Y_size, 1);
        }

        for (arma::uword i = 0; i < Y_size; i++) {
          if constexpr (T == StructType::curves)
          {
            _Y(i, 0) = Rcpp::as<KMA::matrix>(Y[i]);
          }
          else 
          {
            _V(i, 0) = Rcpp::as<KMA::matrix>(Y[i]);
          }
        }
      }

    Rcpp::List probKMA_run()
    {
      bool exe_print = _parameters._exe_print;
      if (exe_print)
        Rcpp::Rcout << "############## STARTING ITERATIONS... ##############" << std::endl;

      /// Set Seed ////////////////////////////////////
      bool set_seed = _parameters._set_seed;
      if (set_seed)
      {
        Rcpp::Environment base("package:base");
        Rcpp::Function setSeed = base["set.seed"];
        unsigned int seed = _parameters._seed;
        setSeed(seed);
      }

      /// Initialization parameters ////////////////////////
      std::size_t iter = 0;
      std::string_view criterion = _parameters._stopCriterion;
      const unsigned int iter_max = _parameters._iter_max;
      const auto transform_function = [this](const KMA::matrix& v0)
        {return std::floor(v0.n_rows * (1 - this->_parameters._max_gap));};
      const double quantile4clean = _parameters._quantile4clean;
      const KMA::vector w = _parameters._w;
      const double alpha = _parameters._alpha;
      const unsigned int iter4clean = _parameters._iter4clean;
      const double tol4clean = _parameters._tol4clean;
      KMA::ivector c = _parameters._c;
      const arma::uword _n_rows_V = _V.n_rows;
      const arma::uword _n_rows_Y = _Y.n_rows;
      Rcpp::Environment stats("package:stats");
      Rcpp::Function quantile = stats["quantile"];

#ifdef _OPENMP
      const unsigned int n_threads = _parameters._n_threads;
#endif

      // Initialization data structures ////////////////////////
      KMA::vector J_iter(iter_max,arma::fill::zeros);
      KMA::vector BC_dist_iter(iter_max,arma::fill::zeros);
      auto BC_dist = std::numeric_limits<double>::infinity();
      KMA::vector sd(2);
      KMA::matrix _D(_n_rows_Y,_n_rows_V);
      std::vector<arma::urowvec> V_dom(_n_rows_V);
      KMA::ivector c_k(_n_rows_V);
      KMA::matrix P_old(_n_rows_Y,_n_rows_V);
      KMA::umatrix D0(_n_rows_Y,_n_rows_V);
      KMA::matrix temp_DP(_n_rows_Y,_n_rows_V);
      KMA::umatrix keep;
      KMA::matrix P(_P0);
      KMA::imatrix S(_S0);

      /// Iterate ////////////////////////////////////
      while(iter < iter_max and BC_dist > _parameters._tol)
      {
        iter++;
        if (exe_print)
          Rcpp::Rcout<<"Iter = "<<iter<<std::endl;

        ///// clean motifs //////////////////////////
        P_old = P;
        if((iter>1)&&(!(iter%iter4clean))&&(BC_dist<tol4clean))
        {
          keep = _D < Rcpp::as<double>(quantile(_D,quantile4clean));
          const KMA::uvector & empty_k = arma::find(arma::sum(keep,0)==0);
          for(arma::uword k : empty_k)
            keep(arma::index_min(_D.col(k)),k) = 1;
          P.zeros();
          P.elem(arma::find(keep==1)).fill(1); // set one where values of keep are true
        }

        ///// compute motifs ///////////////////////
        if((iter != 1) || (!init_motifs)){ 
          for(arma::uword i = 0;i < _n_rows_V;++i)
          {
            const arma::urowvec& V_dom_temp = util::findDomain<KMA::matrix>(_V(i,0));

            const auto& V_new_variant = _motfac->compute_motif(V_dom_temp,
                                                              S.col(i),
                                                              P.col(i),_Y,
                                                              _parameters._m);
            if(auto ptr_1 = std::get_if<MotifPure::indexField>(&V_new_variant))
            {
              Rcpp::Rcout << "Warning: no test for this in compute motifs shifts changes" <<std::endl;
              const arma::sword& index = ptr_1->second;
              S.col(i) += index;
              _V.row(i) = ptr_1->first;
            }
            else
            {
              _V.row(i) = *(std::get_if<KMA::Mfield>(&V_new_variant));
            }
          }
        }

        for(arma::uword i = 0;i < _n_rows_V;++i)
        {
          V_dom[i] = util::findDomain<KMA::matrix>(_V(i,0));
        }

        if((iter>1)&&(!(iter%_parameters._iter4elong))&&(BC_dist<_parameters._tol4elong))
        {
          if (exe_print)
            Rcpp::Rcout<<"Trying to elongate at iter:"<<iter<<std::endl;
          _motfac -> elongate_motifs(_V,V_dom,S,P,
                                     _Y,_D, _parameters,
                                     _perfac,_dissfac,quantile);
        }


        ////// find shift warping minimizing dissimilarities /////////////
        for(arma::uword i = 0;i<_n_rows_V;++i)
        {
          c_k(i) = transform_function(_V(i,0));
          c_k(i) = std::max(c_k(i), c(i));
        }

        #ifdef _OPENMP
          #pragma omp parallel for collapse(2) firstprivate(sd) num_threads(n_threads)
        #endif
        for (arma::uword i = 0; i < _n_rows_V; ++i)
          for (arma::uword j = 0; j < _n_rows_Y; ++j){
            sd = _dissfac->find_diss(_Y.row(j),_V.row(i),w,alpha,c_k(i));
            S(j,i) = sd(0);
            _D(j,i) = sd(1);
          }

        /// compute memberships /////////////////////
        D0 = (_D == 0);
        const KMA::uvector& mult_assign = arma::find(arma::sum(D0,1) > 1);
        for (arma::uword i : mult_assign) {
          Rcpp::warning("Curve has dissimilarity 0 from two different motifs. Using only one of them...");
          const KMA::uvector& indexes = arma::find(D0.row(i) == 1);
          D0.row(i).zeros();
          D0(i,indexes(arma::randi(arma::distr_param(0, indexes.n_elem - 1)))) = 1;
        }

        const KMA::uvector & D0_index = arma::find(arma::sum(D0,1) == 1);
        for(arma::uword i : D0_index) {
          const KMA::uvector & col = arma::find(D0.row(i)==1);
          P(i,col(0)) = 1;
        }

        const KMA::uvector & not_D0_index = arma::find(arma::sum(D0,1) !=1);
        const KMA::matrix & Dm = arma::pow(_D.rows(not_D0_index),1/(_parameters._m-1));
        P.rows(not_D0_index) = 1 / (Dm % arma::repmat(arma::sum(1/Dm,1),1,Dm.n_cols));
        const KMA::uvector & deg_indexes = arma::find(arma::sum(P,0)==0);
        for (arma::uword k : deg_indexes) {
          Rcpp::warning("Motif is degenerate (zero membership). Selecting a new center...");
          P(index_min(_D.col(k)),k) = 1;
        }

        /// evaluate objective functions ////////////////
        temp_DP = _D % (arma::pow(P,_parameters._m));
        temp_DP.replace(arma::datum::nan,0);
        J_iter(iter-1) = arma::accu(temp_DP);

        /// compute Bhattacharyya distance between P_old and P_new ///////////////
        const arma::colvec & BC_dist_k = -arma::log(arma::sum(arma::sqrt(P_old % P),1));
        if (criterion == "max")
          BC_dist = arma::max(BC_dist_k);
        else if (criterion == "mean")
          BC_dist = arma::mean(BC_dist_k);
        else if (criterion == "quantile")
          BC_dist = Rcpp::as<double>(quantile(BC_dist_k,_parameters._quantile));

        BC_dist_iter(iter-1) = BC_dist;
        if (exe_print)
          Rcpp::Rcout<<"BC_dist="<<BC_dist<<std::endl;

      }

      if(iter == iter_max)
        Rcpp::warning("maximum number of iterations reached, method stops");

      /////  prepare output //////////////////////////////////
      _P_clean.set_size(_n_rows_Y,_n_rows_V);
      _P_clean.fill(0);
      _S_clean = S;
      KMA::matrix  D_clean(_n_rows_Y,_n_rows_V);

      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto& pair_motif_shift = _motfac->compute_motif(V_dom[k], S.col(k),
                                                              P.col(k), _Y,
                                                              _parameters._m);
        _V.row(k) = *(std::get_if<KMA::Mfield>(&pair_motif_shift));
      }

      keep = _D < Rcpp::as<double>(quantile(_D,quantile4clean));
      const KMA::uvector& empty_k = arma::find(arma::sum(keep,0) == 0);

      for (arma::uword k: empty_k)
        keep(arma::index_min(_D.col(k)),k) = 1;

      _P_clean(arma::find(keep==1)).fill(1);
      KMA::Mfield V_clean(_n_rows_V,_V.n_cols);
      std::map<arma::sword,arma::sword> shift_s;
      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto& new_motif =  _motfac->compute_motif(V_dom[k], S.col(k),
                                                        arma::conv_to<KMA::vector>::from(_P_clean.col(k)),
                                                         _Y,_parameters._m);
        if (auto ptr = std::get_if<KMA::Mfield>(&new_motif)){
          V_clean.row(k) = *ptr;
        } else {
          Rcpp::Rcout << "Warning: no test for this in compute motifs shifts changes" <<std::endl;
          const auto& pair_motif_shift = std::get_if<std::pair<KMA::Mfield,arma::sword>>(&new_motif);
          V_clean.row(k) = pair_motif_shift->first;
          shift_s.insert(std::make_pair(k, pair_motif_shift->second));
        }
      }
      for(auto it = shift_s.begin();it != shift_s.cend(); ++it){
        _S_clean.col(it->first) += it->second;
      }
      std::vector<arma::urowvec> V_dom_new(_n_rows_V);
      for(arma::uword k=0; k < _n_rows_V ; ++k){
        V_dom_new[k] = util::findDomain<KMA::matrix>(V_clean(k,0));
      }

      /// compute dissimilarities from cleaned motifs, fill D_clean //////////////
      _dissfac -> computeDissimilarityClean(D_clean,_S_clean,V_dom_new,V_clean,_Y);

      /// return output //////////////////////////////////////////////////////
      J_iter.resize(iter);
      BC_dist_iter.resize(iter);
      return toR(V_clean,_P_clean,_S_clean,_D,D_clean,J_iter,BC_dist_iter,iter,P,S);

    }


    void set_parameters(const Rcpp::List& newParameters)
    {
      _parameters = newParameters;

      _dissfac -> set_parameters(_parameters);
    }

    void initial_motifs(const Rcpp::List& V_init,std::string_view diss)
    {
      const Rcpp::List& V0_init = V_init[0];
      const Rcpp::List& V1_init = V_init[1];

      if (diss == "H1") {
        handleCaseH1<StructType::motifs>(V0_init, V1_init);
      } else if (diss == "L2") {
        handleCaseL2<StructType::motifs>(V0_init, V1_init);
      }
    }


    void reinit_motifs(const arma::ivec& c,
                       arma::sword d)
    {
      arma::uword K = c.size();
      _V.set_size(K, _isY0 + _isY1);
      for(arma::uword k=0; k < K; ++k){
        if (_isY0) {
          _V(k,0).set_size(c(k),d);
          _V(k,0).fill(arma::fill::zeros);
        }
        if (_isY1) {
          _V(k,1).set_size(c(k),d);
          _V(k,1).fill(arma::fill::zeros);
        }
      }
    }

    void set_P0(const KMA::matrix& P0)
    {
        _P0 = P0;
    }

    void set_S0(const KMA::imatrix& S0)
    {
        _S0 = S0;
    }

    // return Rcpp::List with all the outputs of probKMA
    Rcpp::List toR(const KMA::Mfield& V_clean,
                   const KMA::imatrix& P_clean,
                   const KMA::imatrix& S_clean,
                   const KMA::matrix& _D,
                   const KMA::matrix& D_clean,
                   const KMA::vector& J_iter,
                   const KMA::vector& BC_dist_iter,
                   std::size_t iter,
                   const KMA::matrix& P,
                   const KMA::imatrix& S) const
    {
        // conv to Rcpp::List of V0,V1,V0_clean,V1_clean
        Rcpp::List V0(_V.n_rows);
        Rcpp::List V1(_V.n_rows);
        Rcpp::List V0_clean(_V.n_rows);
        Rcpp::List V1_clean(_V.n_rows);

        for (arma::uword k = 0; k < _V.n_rows; ++k){
          if(_isY0)
          {
            V0[k] = _V(k,0);
            V0_clean[k] = V_clean(k,0);
          }
          if(_isY1)
          {
            V1[k] = _V(k,1);
            V1_clean[k] = V_clean(k,1);
          }
        }

        if (!_parameters._return_options){
          return Rcpp::List::create(Rcpp::Named("V0") = V0,
                                    Rcpp::Named("V1") = V1,
                                    Rcpp::Named("V0_clean") = V0_clean,
                                    Rcpp::Named("V1_clean") = V1_clean,
                                    Rcpp::Named("P") = P,
                                    Rcpp::Named("P_clean") = P_clean,
                                    Rcpp::Named("S") = S,
                                    Rcpp::Named("S_clean") = S_clean,
                                    Rcpp::Named("D") = _D,
                                    Rcpp::Named("D_clean") = D_clean,
                                    Rcpp::Named("iter") = iter,
                                    Rcpp::Named("J_iter") = J_iter,
                                    Rcpp::Named("BC_dist_iter") = BC_dist_iter);
        } else {
          return pushResult(V0,V1,
                            V0_clean,V1_clean,
                            V_clean,P_clean,
                            S_clean,_D,D_clean,
                            J_iter,BC_dist_iter,iter,P,S);
        }
    }

    //Mandatory since we want to return more than 20 elements in a List
    Rcpp::List pushResult(const Rcpp::List V0,
                          const Rcpp::List V1,
                          const Rcpp::List V0_clean,
                          const Rcpp::List V1_clean,
                          const KMA::Mfield & V_clean,
                          const KMA::imatrix & P_clean,
                          const KMA::imatrix & S_clean,
                          const KMA::matrix & _D,
                          const KMA::matrix & D_clean,
                          const KMA::vector & J_iter,
                          const KMA::vector & BC_dist_iter,
                          std::size_t iter,
                          const KMA::matrix & P,
                          const KMA::imatrix & S) const
    {
      Rcpp::CharacterVector names(32);
      Rcpp::List result(32);
      names[0] = "V0";
      result[0] = V0;
      names[1] = "V1";
      result[1] = V1;
      names[2] = "V0_clean";
      result[2] = V0_clean;
      names[3] = "V1_clean";
      result[3] = V1_clean;
      names[4] = "P0";
      result[4] = _P0;
      names[5] = "P_clean";
      result[5] = P_clean;
      names[6] = "S0";
      result[6] = _S0;
      names[7] = "S_clean";
      result[7] = S_clean;
      names[8] = "D";
      result[8] = _D;
      names[9] = "D_clean";
      result[9] = D_clean;
      names[10] = "iter";
      result[10] = iter;
      names[11] = "J_iter";
      result[11] = J_iter;
      names[12] = "BC_dist_iter";
      result[12] = BC_dist_iter;
      names[13] = "standardize";
      result[13] = _parameters._standardize;
      names[14] = "K";
      result[14] = _parameters._K;
      names[15] = "c";
      result[15] = _parameters._c;
      names[16] = "c_max";
      result[16] = _parameters._c_max;
      names[17] = "iter_max";
      result[17] = _parameters._iter_max;
      names[18] = "quantile";
      result[18] = _parameters._quantile;
      names[19] = "tol";
      result[19] = _parameters._tol;
      names[20] = "stop_criterion";
      result[20] = _parameters._stopCriterion;
      names[21] = "m,";
      result[21] = _parameters._m;
      names[22] = "iter4elong";
      result[22] = _parameters._iter4elong;
      names[23] = "tol4elong";
      result[23] = _parameters._tol4elong;
      names[24] = "max_elong";
      result[24] = _parameters._max_elong;
      names[25] = "trials_elong";
      result[25] = _parameters._trials_elong;
      names[26] = "deltaJk_elong";
      result[26] = _parameters._deltaJK_elong;
      names[27] = "max_gap";
      result[27] = _parameters._max_gap;
      names[28] = "iter4clean";
      result[28] = _parameters._iter4clean;
      names[29] = "tol4clean";
      result[29] = _parameters._tol4clean;
      names[30] = "P";
      result[30] = P;
      names[31] = "S";
      result[31] = S;

      result.attr("names") = Rcpp::wrap(names);
      return result;
    }


    Rcpp::List compute_silhouette(bool align)
    { 
      const double alpha = _parameters._alpha;
      const KMA::vector & w = _parameters._w;
      const KMA::matrix & first_y0 = _Y(0,0);
      const arma::uword d = first_y0.n_cols;
      const arma::uword K = _parameters._K;

      std::vector<KMA::uvector> V_dom(K); 
      KMA::ivector V_length(K);
      KMA::uvector v_dom_k;
      for(arma::uword k = 0; k < K; ++k)
      {
        V_length(k) = _V(k,0).n_rows;
        v_dom_k.set_size(_V(k,0).n_rows);
        for (arma::uword j = 0; j < _V(k,0).n_rows; ++j)
        {
          const arma::uvec & nan_v_row_j = arma::find_nan(_V(k,0).row(j));
          v_dom_k[j] = (nan_v_row_j.n_elem != d); 
        }
        V_dom[k] = v_dom_k;
      }
      
      std::vector<KMA::uvector> curves_in_motifs(_P_clean.n_cols); 
      for (arma::uword k = 0; k < K; ++k) 
      {
        const KMA::ivector & P_clean_k = _P_clean.col(k);
        const KMA::uvector & P_clean_1 = find(P_clean_k == 1); 
        curves_in_motifs[k] = P_clean_1;
      }
 
      KMA::ivector curves_in_motifs_number = arma::sum(_P_clean,0).t();

      // for each centroid take the shift of the curve associated to that centroid 
      std::vector<KMA::ivector> S_clean_k(K);
      for (unsigned int k=0; k < K; ++k)
      {
        const KMA::ivector & col_S_clean = _S_clean.col(k);
        S_clean_k[k] = col_S_clean.elem(curves_in_motifs[k]);
      }
 
      // compute distances between pieces of curves
      const arma::uword Y_in_motifs_size = arma::accu(_P_clean);
      KMA::Mfield Y_in_motifs(Y_in_motifs_size,_isY0 + _isY1);
      KMA::uvector index;
      arma::uword l = 0;
      unsigned int y_len;
      for (unsigned int k= 0; k < K; ++k) // for each centroid
      { 
        const KMA::uvector & curves_in_motif = curves_in_motifs[k];
        const KMA::uvector & v_dom = V_dom[k];
        for(unsigned int j = 0; j < curves_in_motif.n_elem; ++j) // for each curve assigned to that centroid
        { 
          const int s = S_clean_k[k][j];
          index = arma::regspace<arma::uvec>(1,v_dom.n_elem - std::max(0,1-s)) + std::max(1,s) - 1;
          int index_size = index.n_elem;
          const KMA::matrix & y0 = _Y(curves_in_motif[j],0);
          y_len = y0.n_rows;
          Y_in_motifs(l,0).set_size(v_dom.n_elem,d);
          Y_in_motifs(l,0).fill(arma::datum::nan);
          auto filtered_j = std::views::iota(0,index_size)  
            | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
          for(int j : filtered_j) 
            Y_in_motifs(l,0).row(std::max(0,1-s) + j) =  y0.row(index[j] - 1);
          if (_isY0 and _isY1)
          {
            const KMA::matrix & y1 = _Y(curves_in_motif[j],1);
            y_len = y1.n_rows;
            Y_in_motifs(l,1).set_size(v_dom.n_elem,d);
            Y_in_motifs(l,1).fill(arma::datum::nan);
            auto filtered_j = std::views::iota(0,index_size) 
              | std::views::filter([&y_len,&v_dom, &index](int j){return (index[j] <= y_len && v_dom(j));});
            for(int j : filtered_j)
              Y_in_motifs(l,1).row(std::max(0,1-s) + j) =  y1.row(index[j] - 1);
          }
          l++;
        }
      }

      
      KMA::uvector Y_motifs = util::repLem<arma::uvec>(arma::regspace<arma::uvec>(0,K-1),curves_in_motifs_number);

      const arma::uword Y_motifs_size = Y_motifs.size();

      KMA::umatrix indeces_YY = util::combn2<arma::uword>(arma::regspace<arma::uvec>(0,Y_in_motifs_size-1)); //combn of indeces of Y_in_motifs

      KMA::ivector V_length_Y_motifs = util::repLem<arma::ivec>(V_length,curves_in_motifs_number);
      
      const KMA::ivector c = _parameters._c;
      
      KMA::ivector c_Y_motifs = util::repLem<arma::ivec>(c,curves_in_motifs_number);
      
      KMA::imatrix YY_lengths = util::combn2<arma::sword>(V_length_Y_motifs);
      
      const int YY_length_size = YY_lengths.n_cols;

      const arma::irowvec & YY_length_row0 = YY_lengths.row(0);

      const arma::irowvec & YY_length_row1 = YY_lengths.row(1);

      const arma::urowvec & swap = (YY_length_row0 < YY_length_row1);

      const arma::urowvec & equal_length = (YY_length_row0 == YY_length_row1);

      auto filtered_j_swap = std::views::iota(0,YY_length_size) 
           | std::views::filter([&swap](int j){return swap(j);});
      
      for (int j : filtered_j_swap)
        std::swap(indeces_YY(0,j),indeces_YY(1,j));
      
      KMA::vector SD(YY_length_size);
      KMA::vector min_diss;
      _dissfac->set_both(true); // utile solo se _transformed = true
      if(align)
      {
        #ifdef _OPENMP
            #pragma omp parallel for firstprivate(min_diss)
        #endif
        for(int i = 0; i < YY_length_size; ++i)
        {
          min_diss = _dissfac->find_diss_aligned(Y_in_motifs.row(indeces_YY(0,i)),    
                                                 Y_in_motifs.row(indeces_YY(1,i)),
                                                 equal_length(i));
                                                                  
          SD(i) = min_diss(1);
        }
      }
      else
      { 
        KMA::imatrix c_Y_motifs_comb = util::combn2<arma::sword>(c_Y_motifs);  
        #ifdef _OPENMP
            #pragma omp parallel for firstprivate(min_diss)
        #endif
        for(int i = 0; i < YY_length_size; ++i){
          min_diss = _dissfac->find_diss(Y_in_motifs.row(indeces_YY(0,i)),
                                         Y_in_motifs.row(indeces_YY(1,i)),
                                         w, alpha, std::min(c_Y_motifs_comb(0,i),
                                                            c_Y_motifs_comb(1,i)));
                                     
          SD(i) = min_diss(1);
        }
      }

      KMA::matrix YY_D(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
      arma::uword k = 0;
      for (arma::uword j = 0; j < Y_motifs_size; ++j){
        for (arma::uword i = j+1; i < Y_motifs_size; ++i){
          YY_D(i,j) = SD(k);
          k++;
        }
      }
      
      YY_D = YY_D + YY_D.t(); 
        
      // compute intra-cluster distances
      KMA::umatrix intra(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
      for(unsigned int i = 0; i < Y_motifs_size; ++i){
        const KMA::uvector & temp = (Y_motifs == Y_motifs(i));
        intra.col(i) = temp;
      }
      intra.diag().fill(0);
      
      const KMA::vector & curves_in_motifs_number_Y_motifs = arma::conv_to<KMA::vector>::from(curves_in_motifs_number.elem(Y_motifs));
      const KMA::vector & a = arma::sum(intra%YY_D, 0).t()/(curves_in_motifs_number_Y_motifs - 1);

      // compute inter-cluster distances
      Y_motifs += 1;
      KMA::uvector Y_motifs_mod(Y_motifs_size);
      for(unsigned int i=0; i < Y_motifs_size; ++i)
          Y_motifs_mod(i) = Y_motifs(i)+1>K ? (Y_motifs(i)+1)%K - 1 : Y_motifs(i);
      const KMA::vector & curves_in_motifs_number_rep = arma::conv_to<KMA::vector>::from(curves_in_motifs_number.elem(Y_motifs_mod));
      
      KMA::umatrix inter(Y_motifs_size,Y_motifs_size,arma::fill::zeros);
      KMA::matrix b_k(K-1,Y_motifs_size,arma::fill::zeros);
      arma::uword motif;
      arma::uword inter_motif;
      for(unsigned int k = 1; k <= K-1; ++k){
        for(unsigned int i = 0; i < Y_motifs_size; ++i){
          motif = Y_motifs(i);
          inter_motif = (motif+k)>K? (motif+k)%K : motif+k; 
          inter.col(i) = (Y_motifs== inter_motif);
        }
        b_k.row(k-1) = sum(inter%YY_D,0)/(curves_in_motifs_number_rep.t());
      }
      
      arma::rowvec b_tmp(Y_motifs_size);
      if(b_k.n_rows > 1){
        b_tmp = arma::min(b_k,0);
      }else{
        b_tmp = b_k;
      }
      arma::vec b = b_tmp.t();
      
      // compute silhouette
      KMA::vector silhouette= (b-a)/arma::max(a,b);
      for(auto & sil : silhouette){ 
        if(!arma::is_finite(sil))
          sil = 0;
      }
      
      // compute average silhouette per cluster
      KMA::vector silhouette_average(K);
      silhouette_average.fill(arma::datum::nan);
      KMA::vector silhouette_k;
      KMA::uvector indexes;
      KMA::uvector sorted_indexes;
      for (unsigned int k = 0; k < K; k++) { 
        indexes = find(Y_motifs == k + 1);
        silhouette_k = silhouette.elem(indexes);
        sorted_indexes = arma::sort_index(silhouette_k, "descend");
        curves_in_motifs[k] = curves_in_motifs[k].elem(sorted_indexes);
        curves_in_motifs[k] += 1;
        silhouette.elem(indexes) = arma::sort(silhouette_k, "descend");
        silhouette_average(k) = mean(silhouette_k);
      }

      _dissfac->set_both(false); // per usi futuri 

      return Rcpp::List::create(Rcpp::Named("silhouette") = silhouette,
                                Rcpp::Named("motifs") = Y_motifs,
                                Rcpp::Named("curves") = curves_in_motifs,
                                Rcpp::Named("silhouette_average") = silhouette_average,
                                Rcpp::Named("curves_in_motifs_number") = curves_in_motifs_number);

    }

    // Functional Data
    KMA::Mfield _Y;
    KMA::Mfield _V;

    //Motif and dissimilarity
    std::shared_ptr<MotifPure> _motfac;
    std::shared_ptr<Dissimilarity> _dissfac;
    std::shared_ptr<PerformanceIndexAB> _perfac;

    //Parameters
    Parameters _parameters;

    // Membership and shifting matrix
    KMA::matrix _P0;
    KMA::imatrix _S0;
    KMA::imatrix _P_clean;
    KMA::imatrix _S_clean;
    bool _isY0 = true;
    bool _isY1 = true;
    bool init_motifs = false;
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,
                 const Rcpp::List& parameters,
                 const KMA::matrix& P0,const KMA::imatrix& S0,
                 const std::string& diss):
                 _probKMA(std::make_unique<_probKMAImp>(Y,parameters,P0,S0,diss)) {};

ProbKMA::ProbKMA(const Rcpp::List& Y,
                 const Rcpp::List& parameters,
                 const KMA::matrix& P0,const KMA::imatrix& S0,
                 const std::string& diss,
                 const Rcpp::List & V_init):
                 _probKMA(std::make_unique<_probKMAImp>(Y,parameters,P0,S0,diss,V_init)) {};


Rcpp::List ProbKMA::probKMA_run() const
{
   return _probKMA -> probKMA_run();
}

void ProbKMA::set_parameters(const Rcpp::List& newParameters)
{
    _probKMA -> set_parameters(newParameters);
}

void ProbKMA::reinit_motifs(const arma::ivec & c,
                            arma::sword d)
{
    _probKMA -> reinit_motifs(c,d);
}

void ProbKMA::set_P0(const KMA::matrix& P0)
{
    _probKMA -> set_P0(P0);
}

void ProbKMA::set_S0(const KMA::imatrix& S0)
{
    _probKMA -> set_S0(S0);
}

Rcpp::List ProbKMA::compute_silhouette(bool align)
{
  return _probKMA -> compute_silhouette(align);
}


RCPP_EXPOSED_CLASS(ProbKMA);

RCPP_MODULE(ProbKMAModule) {
  Rcpp::class_<ProbKMA>("ProbKMA")
  .constructor<Rcpp::List,Rcpp::List,
               KMA::matrix,KMA::imatrix,std::string>()
  .constructor<Rcpp::List,Rcpp::List,
               KMA::matrix,KMA::imatrix,
               std::string,Rcpp::List>()
  .method("probKMA_run",&ProbKMA::probKMA_run)
  .method("set_parameters", &ProbKMA::set_parameters)
  .method("reinit_motifs", &ProbKMA::reinit_motifs)
  .method("set_P0", &ProbKMA::set_P0)
  .method("set_S0", &ProbKMA::set_S0)
  .method("compute_silhouette", &ProbKMA::compute_silhouette);
}
