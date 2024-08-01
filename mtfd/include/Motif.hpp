#ifndef MOTIF_HPP
#define MOTIF_HPP
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include <Rcpp.h>
#include <numeric>
#include <ranges>
#include <algorithm>
#include <variant>
#include <memory>
#include "Utilities.hpp"
#include "Parameters.hpp"
#include "PerformanceIndex.hpp"
#include "Dissimilarity.hpp"


class MotifPure
{
  public:
    using indexField = std::pair<KMA::Mfield,arma::sword>;
    
    MotifPure() = default;

    virtual std::variant<indexField,KMA::Mfield>
    compute_motif(const arma::urowvec& v_dom,
                  const KMA::ivector& s_k,
                  const KMA::vector& p_k,
                  const KMA::Mfield& Y,
                  double m) const = 0;
    
    virtual void set_parameters(const Parameters & newParameters) = 0;

    virtual
    void elongate_motifs(KMA::Mfield& V_new,
                         std::vector<arma::urowvec>& V_dom,
                         KMA::imatrix& S_k,const KMA::matrix& P_k,
                         const KMA::Mfield& Y,const KMA::matrix& D,
                         const Parameters& param,
                         const std::shared_ptr<PerformanceIndexAB>& perf,
                         const std::shared_ptr<Dissimilarity>& diss,
                         const Rcpp::Function & quantile_func) const = 0;
    
                          
    virtual ~MotifPure() = default;
    
};

class MotifSobol: public MotifPure
{
public:
  
  MotifSobol(bool transformed);

  virtual ~MotifSobol() override = default;
  
protected:
  
  KMA::matrix compute_v_new(const KMA::Mfield& Y_inters_k,
                            const KMA::umatrix& Y_inters_supp,
                            const arma::urowvec & v_dom,
                            arma::uword v_len,
                            const KMA::vector & p_k,
                            arma::uword d,
                            arma::uword m) const;
  // Use1 = True --> H1
  template<bool use1>
  std::variant<indexField,KMA::Mfield> 
  compute_motif_helper(const arma::urowvec& v_dom,
                       const KMA::ivector& s_k,
                       const KMA::vector& p_k,
                       const KMA::Mfield& Y,
                       double m) const;
  
  template<bool use1>
  void elongate_motifs_helper(KMA::Mfield& V_new,
                              std::vector<arma::urowvec>& V_dom,
                              KMA::imatrix& S_k,const KMA::matrix& P_k,
                              const KMA::Mfield& Y,const KMA::matrix& D,
                              const Parameters& param,
                              const std::shared_ptr<PerformanceIndexAB>& perf,
                              const std::shared_ptr<Dissimilarity>& diss,
                              const Rcpp::Function & quantile_func) const;
  
  template<bool use1>
  void elongation(KMA::Mfield& V_new, 
                  std::vector<arma::urowvec> & V_dom,  
                  KMA::imatrix & S_k, 
                  const arma::vec & p_k, 
                  const arma::ivec& len_elong_k, 
                  const arma::uvec& keep_k,
                  double c,
                  const KMA::Mfield Y, 
                  const unsigned int index,
                  const Parameters& param,
                  const std::shared_ptr<PerformanceIndexAB>& performance,
                  const std::shared_ptr<Dissimilarity>& diss) const;

  bool _transformed;
};


class MotifL2 final: public MotifSobol
{
public:
  
  MotifL2(bool transformed);
  
  virtual std::variant<indexField,KMA::Mfield>
    compute_motif(const arma::urowvec& v_dom,
                  const KMA::ivector& s_k,
                  const KMA::vector& p_k,
                  const KMA::Mfield& Y,
                  double m) const override;

  void set_parameters(const Parameters & newParameters) override;
  
  virtual
  void elongate_motifs(KMA::Mfield& V_new,
                       std::vector<arma::urowvec>& V_dom,
                       KMA::imatrix& S_k,const KMA::matrix& P_k,
                       const KMA::Mfield& Y,const KMA::matrix& D,
                       const Parameters& param,
                       const std::shared_ptr<PerformanceIndexAB>& perf,
                       const std::shared_ptr<Dissimilarity>& diss,
                       const Rcpp::Function & quantile_func) const override;
  
  virtual ~MotifL2() override = default;
  
};

class MotifH1 final: public MotifSobol
{
public:

  MotifH1(bool transformed);
  
  virtual ~MotifH1() override = default;
  
  virtual std::variant<indexField,KMA::Mfield>
    compute_motif(const arma::urowvec& v_dom,
                  const KMA::ivector& s_k,
                  const KMA::vector& p_k,
                  const KMA::Mfield& Y,
                  double m) const override;

  void set_parameters(const Parameters & newParameters) override;
  
  virtual
  void elongate_motifs(KMA::Mfield& V_new,
                       std::vector<arma::urowvec>& V_dom,
                       KMA::imatrix& S_k,const KMA::matrix& P_k,
                       const KMA::Mfield& Y,const KMA::matrix& D,
                       const Parameters& param,
                       const std::shared_ptr<PerformanceIndexAB>& perf,
                       const std::shared_ptr<Dissimilarity>& diss,
                       const Rcpp::Function & quantile_func) const override;
  
};

#include "Motif.ipp"

#endif // MOTIF_HPP

