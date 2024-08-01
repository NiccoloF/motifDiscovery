#ifndef DISSIMILARITY_HPP
#define DISSIMILARITY_HPP
#include "TypeTraits.hpp"
#include "Parameters.hpp"
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <Utilities.hpp>
#include <ranges>
#include "Parameters.hpp"

// Abstract class for dissimilarities
class Dissimilarity
{
public:

  Dissimilarity() = default;

  virtual ~Dissimilarity() = default;

  virtual void set_both(bool both) = 0;

  // compute dissimilarity
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const = 0;

  virtual void set_parameters(const Parameters & newParameters) = 0;

  // compute dissimilarity from cleaned motifs
  virtual void computeDissimilarityClean(KMA::matrix & D_clean,
                                         const KMA::imatrix & S_clean,
                                         const std::vector<arma::urowvec> & V_dom_new,
                                         const KMA::Mfield & V_clean,
                                         const KMA::Mfield & Y) const = 0;

  // Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w,
                                double alpha, unsigned int c_k) const = 0;

  virtual KMA::vector find_diss_aligned(const KMA::Mfield Y,
                                        const KMA::Mfield V,
                                        bool aligned) const = 0;

protected:

  virtual double distance(const KMA::matrix& y,
                          const KMA::matrix& v) const = 0;
};

class SobolDiss : public Dissimilarity
{
public:

    SobolDiss(const KMA::vector& w,
              bool transformed);

protected:

    virtual double distance(const KMA::matrix& y,
                            const KMA::matrix& v) const override;

    template<bool use1>
    KMA::vector find_diss_helper(const KMA::Mfield Y,
                                 const KMA::Mfield V,
                                 const KMA::vector& w,
                                 double alpha, unsigned int c_k) const;

    template<bool use1>
    KMA::vector find_diss_aligned_helper(const KMA::Mfield Y,
                                         const KMA::Mfield V,
                                         bool aligned) const;

    template<bool use1>
    void computeDissimilarityClean_helper(KMA::matrix & D_clean,
                                          const KMA::imatrix & S_clean,
                                          const std::vector<arma::urowvec> & V_dom_new,
				                                  const KMA::Mfield & V_clean,
                                          const KMA::Mfield & Y) const;

    KMA::vector _w;
    bool _transformed;
    bool _both = false;

};

#include "Dissimilarity.ipp"

class L2 final: public SobolDiss
{
public:

  L2(const KMA::vector& w, bool transformed);
  virtual ~L2() = default;

  virtual void set_both(bool both) override;

  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;

  virtual void computeDissimilarityClean(KMA::matrix & D_clean,
                                         const KMA::imatrix & S_clean,
                                         const std::vector<arma::urowvec> & V_dom_new,
                                         const KMA::Mfield & V_clean,
                                         const KMA::Mfield & Y) const override;

  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w,
                                double alpha, unsigned int c_k) const override;

  virtual KMA::vector find_diss_aligned(const KMA::Mfield Y,
                                        const KMA::Mfield V,
                                        bool aligned) const override;

  void set_parameters(const Parameters & newParameters) override;

};

class H1 final: public SobolDiss
{
public:

  H1(const KMA::vector& w,double alpha, bool transformed);
  virtual ~H1() = default;

  virtual void set_both(bool both) override;

  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;

  void set_parameters(const Parameters & newParameters) override;

  virtual void computeDissimilarityClean(KMA::matrix & D_clean,
                                         const KMA::imatrix & S_clean,
                                         const std::vector<arma::urowvec> & V_dom_new,
                                         const KMA::Mfield & V_clean,
                                         const KMA::Mfield & Y) const override;

  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w,
                                double alpha, unsigned int c_k) const override;

  virtual KMA::vector find_diss_aligned(const KMA::Mfield Y,
                                        const KMA::Mfield V,
                                        bool aligned) const override;

  double _alpha;

};




#endif // DISSIMILARITY_HPP
