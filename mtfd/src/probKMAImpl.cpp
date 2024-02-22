#include "ProbKMA.hpp"
#include "Factory.hpp"
#include <forward_list>
#include <limits>
#include <string_view>

class ProbKMA::_probKMAImp
{
public:

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
                    dissfac.FactoryRegister<L2>("L2",_parameters._w);
                    dissfac.FactoryRegister<H1>("H1",_parameters._w,_parameters._alpha);

                    //Create Motif factory
                    util::SharedFactory<MotifPure> motfac;
                    motfac.FactoryRegister<MotifL2>("L2");
                    motfac.FactoryRegister<MotifH1>("H1");

                    //Create Performance factory
                    util::SharedFactory<PerformanceIndexAB> perfac;
                    perfac.FactoryRegister<PerformanceL2>("L2");
                    perfac.FactoryRegister<PerformanceH1>("H1");

                    //check whether diss is valid
                    _motfac = motfac.instantiate(diss); // copy elision
                    _dissfac = dissfac.instantiate(diss); // copy elision
                    _perfac = perfac.instantiate(diss); // copy elision
                    if(not(_motfac and _dissfac and _perfac))
                      Rcpp::Rcerr<<"Invalid dissimilarity: Choose between L2,H1"<<std::endl;
                }

    ~_probKMAImp() = default;

    void Initialize(const Rcpp::List& Y,std::string_view diss)
    {
      // Convert Rcpp Data Structure(List) into Armadillo data Structure(field)
      const Rcpp::List& Y0 = Y[0];
      const Rcpp::List& Y1 = Y[1];

      if (diss == "H1") {
        handleCaseH1(Y0, Y1);
      } else if (diss == "L2") {
        handleCaseL2(Y0, Y1);
      }
      reinit_motifs(_parameters._c,_Y.front().n_cols);
    }

      // Support function for the case "H1"
      void handleCaseH1(const Rcpp::List& Y0, const Rcpp::List& Y1) {
        const arma::uword Y_size = Y0.size();
        _Y.set_size(Y_size, 2);

        for (arma::uword i = 0; i < Y_size; i++) {
          _Y(i, 0) = Rcpp::as<KMA::matrix>(Y0[i]);
          _Y(i, 1) = Rcpp::as<KMA::matrix>(Y1[i]);
        }
      }

      // Support function for the case "L2"
      void handleCaseL2(const Rcpp::List& Y0, const Rcpp::List& Y1) {
        if (!Rf_isNull(Y0[0])) {
          _isY1 = false;
          handleNonNullY(Y0);
        } else {
          _isY0 = false;
          handleNonNullY(Y1);
        }
      }

      // Support function for the case "L2"
      void handleNonNullY(const Rcpp::List& Y) {
        const arma::uword Y_size = Y.size();
        _Y.set_size(Y_size, 1);

        for (arma::uword i = 0; i < Y_size; i++) {
          _Y(i, 0) = Rcpp::as<KMA::matrix>(Y[i]);
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
        for(arma::uword i = 0;i < _n_rows_V;++i)
        {
          const arma::urowvec& V_dom_temp = util::findDomain<KMA::matrix>(_V(i,0));

          const auto& V_new_variant = _motfac->compute_motif(V_dom_temp,
                                                             S.col(i),
                                                             P.col(i),_Y,
                                                             _parameters._m);
          if(auto ptr_1 = std::get_if<MotifPure::indexField>(&V_new_variant))
          {
            const arma::sword& index = ptr_1->second;
            S.col(i) += index;
            _V.row(i) = ptr_1->first;
          }
          else
          {
            _V.row(i) = *(std::get_if<KMA::Mfield>(&V_new_variant));
          }
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
      KMA::matrix  P_clean(_n_rows_Y,_n_rows_V,arma::fill::zeros);
      KMA::imatrix S_clean(S);
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

      P_clean(arma::find(keep==1)).fill(1);
      KMA::Mfield V_clean(_n_rows_V,_V.n_cols);
      std::map<arma::sword,arma::sword> shift_s;
      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto& new_motif =  _motfac->compute_motif(V_dom[k], S.col(k),
                                                         P_clean.col(k), _Y,
                                                         _parameters._m);
        if (auto ptr = std::get_if<KMA::Mfield>(&new_motif)){
          V_clean.row(k) = *ptr;
        } else {
          const auto& pair_motif_shift = std::get_if<std::pair<KMA::Mfield,arma::sword>>(&new_motif);
          V_clean.row(k) = pair_motif_shift->first;
          shift_s.insert(std::make_pair(k, pair_motif_shift->second));
        }
      }
      for(auto it = shift_s.begin();it != shift_s.cend(); ++it){
        S_clean.col(it->first) += it->second;
      }
      std::vector<arma::urowvec> V_dom_new(_n_rows_V);
      for(arma::uword k=0; k < _n_rows_V ; ++k){
        V_dom_new[k] = util::findDomain<KMA::matrix>(V_clean(k,0));
      }

      /// compute dissimilarities from cleaned motifs, fill D_clean //////////////
      _dissfac -> computeDissimilarityClean(D_clean,S_clean,V_dom_new,V_clean,_Y);

      /// return output //////////////////////////////////////////////////////
      J_iter.resize(iter);
      BC_dist_iter.resize(iter);
      return toR(V_clean,P_clean,S_clean,_D,D_clean,J_iter,BC_dist_iter,iter,P,S);

    }


    void set_parameters(const Rcpp::List& newParameters)
    {
      _parameters = newParameters;

      _dissfac -> set_parameters(_parameters);
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

    Rcpp::List get_motifs() const{
      Rcpp::List V0(_V.n_rows);
      Rcpp::List V1(_V.n_rows);
      for(arma::uword k = 0; k < _V.n_rows; ++k)
      {
        if (_isY0)
         V0[k] = _V(k,0);
        if (_isY1)
         V1[k] = _V(k,1);
      }
      return Rcpp::List::create(V0,V1);
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
                   const KMA::matrix& P_clean,
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
                          const KMA::matrix & P_clean,
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
    bool _isY0 = true;
    bool _isY1 = true;
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,
                 const Rcpp::List& parameters,
                 const KMA::matrix& P0,const KMA::imatrix& S0,
                 const std::string& diss):
                 _probKMA(std::make_unique<_probKMAImp>(Y,parameters,P0,S0,diss)) {};


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


Rcpp::List ProbKMA::get_motifs() const
{
  return _probKMA -> get_motifs();
}

// [[Rcpp::export(initialChecks)]]
Rcpp::List initialChecks(const Rcpp::List& Y0,const Rcpp::List& Y1,
                         const Rcpp::NumericMatrix& P0,
                         const Rcpp::NumericMatrix& S0,
                         const Rcpp::List& params,
                         const Rcpp::String& diss,
                         long long seed)
{
  try
  {
    Rcpp::Function checks = Rcpp::Environment::namespace_env("ProbKMAcpp")[".initialChecks"];
    // Call R checks and updata data and parameters
    return checks(Y0,Y1,P0,S0,params,diss,seed);
  }catch (Rcpp::exception& e) {
    // Handle the Rcpp exception
    Rcpp::Rcerr << "Caught exception: " << e.what() << std::endl;
    return Rcpp::List::create();

  } catch (...){
    // Handle other types of exceptions
    Rcpp::Rcerr << "Caught unknown exception." << std::endl;
    return Rcpp::List::create();
  }
}

RCPP_EXPOSED_CLASS(ProbKMA);

RCPP_MODULE(ProbKMAModule) {
  Rcpp::class_<ProbKMA>("ProbKMA")
  .constructor<Rcpp::List,Rcpp::List,
               KMA::matrix,KMA::imatrix,std::string>()
  .method("probKMA_run",&ProbKMA::probKMA_run)
  .method("set_parameters", &ProbKMA::set_parameters)
  .method("get_motifs", &ProbKMA::get_motifs)
  .method("reinit_motifs", &ProbKMA::reinit_motifs)
  .method("set_P0", &ProbKMA::set_P0)
  .method("set_S0", &ProbKMA::set_S0);
}
