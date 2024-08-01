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

    ProbKMA(const Rcpp::List& Y,
            const Rcpp::List& parameters,
            const KMA::matrix& P0,const KMA::imatrix& S0,
            const std::string& diss, 
            const Rcpp::List& V_init);

    // Y: a list containing two list -> Y0 and Y1
    ProbKMA(const Rcpp::List& Y,
            const Rcpp::List& parameters,
            const KMA::matrix& P0,const KMA::imatrix& S0,
            const std::string& diss);

    virtual ~ProbKMA() = default;

    // run probKMA's algorithm
    Rcpp::List probKMA_run() const;

    void set_parameters(const Rcpp::List& parameters);

    void reinit_motifs(const arma::ivec& c, arma::sword d);

    void set_P0(const KMA::matrix& P0);

    void set_S0(const KMA::imatrix& S0);

    Rcpp::List compute_silhouette(bool align); 
    
 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<_probKMAImp> _probKMA;
};

 
#endif // PROBKMA_HPP
