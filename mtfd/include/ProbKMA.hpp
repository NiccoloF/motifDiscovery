#ifndef PROBKMA_HPP
#define PROBKMA_HPP

#include "RcppArmadillo.h"
#include "Parameters.hpp"
#include "Dissimilarity.hpp"
#include "Motif.hpp"
#include "TypeTraits.hpp"
#include <vector>
#include <ranges>
#include <algorithm>
#include <memory>
#include <Rcpp.h>


// Forward declaration
class _probKMAImp;

class ProbKMA
{
 public:

    // Y: a list containing two list -> Y0 and Y1
    ProbKMA(const Rcpp::List& Y,
            const Rcpp::List& parameters,
            const KMA::matrix& P0,const KMA::imatrix& S0,
            const std::string& diss);

    virtual ~ProbKMA() = default;

    // run probKMA's algorithm
    Rcpp::List probKMA_run() const;

    void set_parameters(const Rcpp::List& parameters);

    Rcpp::List get_motifs() const;

    void reinit_motifs(const arma::ivec& c, arma::sword d);

    void set_P0(const KMA::matrix& P0);

    void set_S0(const KMA::imatrix& S0);

 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<_probKMAImp> _probKMA;
};

  Rcpp::List initialChecks(const Rcpp::List& Y0,const Rcpp::List& Y1,
                           const Rcpp::NumericMatrix& P0,
                           const Rcpp::NumericMatrix& S0,
                           const Rcpp::List& params,
                           const Rcpp::String diss);


#endif // PROBKMA_HPP
