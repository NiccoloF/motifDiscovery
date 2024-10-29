#include "Parameters.hpp"

 Parameters::Parameters(const Rcpp::List& params)
  {
    _standardize = params["standardize"];
    _K = params["K"];
    _c = Rcpp::as<arma::ivec>(params["c"]);
    _c_max = Rcpp::as<arma::ivec>(params["c_max"]);
    _iter_max = params["iter_max"];
    _quantile = params["quantile"];
    _stopCriterion = Rcpp::as<std::string>(params["stopCriterion"]);
    _m = params["m"];
    _w = Rcpp::as<arma::vec>(params["w"]);
    _alpha = params["alpha"];
    _tol = params["tol"];
    _iter4elong = params["iter4elong"];
    _tol4elong = params["tol4elong"];
    _max_elong = params["max_elong"];
    _trials_elong = params["trials_elong"];
    _deltaJK_elong = params["deltaJK_elong"];
    _max_gap = params["max_gap"];
    _iter4clean = params["iter4clean"];
    _tol4clean = params["tol4clean"];
    _quantile4clean = params["quantile4clean"];
    _return_options = params["return_options"];
    _seed = params["seed"];
    _exe_print = params["exe_print"];
    _set_seed = params["set_seed"];
    _n_threads = params["n_threads"];
    _transformed = params["transformed"];
  }

Rcpp::List Parameters::to_list()
{
  Rcpp::CharacterVector names(26);
  Rcpp::List result(26);
  names[0] = "standardize";
  result[0] = _standardize;
  names[1] = "K";
  result[1] = _K;
  names[2] = "c";
  result[2] = _c;
  names[3] = "c_max";
  result[3] = _c_max;
  names[4] = "iter_max";
  result[4] = _iter_max;
  names[5] = "quantile";
  result[5] = _quantile;
  names[6] = "stopCriterion";
  result[6] = _stopCriterion;
  names[7] = "m";
  result[7] = _m;
  names[8] = "w";
  result[8] = _w;
  names[9] = "alpha";
  result[9] = _alpha;
  names[10] = "tol";
  result[10] = _tol;
  names[11] = "iter4elong";
  result[11] = _iter4elong;
  names[12] = "tol4elong";
  result[12] = _tol4elong;
  names[13] = "max_elong";
  result[13] = _max_elong;
  names[14] = "trials_elong";
  result[14] = _trials_elong;
  names[15] = "return_options";
  result[15] = _max_elong;
  names[16] = "quantile4clean";
  result[16] = _trials_elong;
  names[17] = "deltaJk_elong";
  result[17] = _deltaJK_elong;
  names[18] = "max_gap";
  result[18] = _max_gap;
  names[19] = "iter4clean";
  result[19] = _iter4clean;
  names[20] = "tol4clean";
  result[20] = _tol4clean;
  names[21] = "seed";
  result[21] = _seed;
  names[22] = "exe_print";
  result[22] = _exe_print;
  names[23] = "set_seed";
  result[23] = _set_seed;
  names[24] = "n_threads";
  result[24] = _n_threads;
  names[25] = "transformed";
  result[25] = _transformed;
  result.attr("names") = Rcpp::wrap(names);
  return result;
}
