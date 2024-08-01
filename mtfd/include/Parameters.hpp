#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <string>
#include "RcppArmadillo.h"


struct Parameters
{
    Parameters() = delete;
    Parameters(const Parameters&) = default;
    Parameters(const Rcpp::List& params);
    Rcpp::List to_list();

    bool _standardize;
    unsigned int _K;
    arma::ivec _c;
    arma::ivec _c_max;
    unsigned int _iter_max;
    double _quantile;
    std::string _stopCriterion;
    double _m;
    arma::vec _w;
    double _alpha;
    double _tol;
    unsigned int _iter4elong;
    double _tol4elong;
    double _max_elong;
    unsigned int _trials_elong;
    double _deltaJK_elong;
    double _max_gap;
    unsigned int _iter4clean;
    double _tol4clean;
    double _quantile4clean;
    bool _return_options;
    unsigned int _seed;
    bool _exe_print;
    bool _set_seed;
    unsigned int _n_threads;
    bool _transformed;
};

#endif // PARAMETERS_HPP
