#include "Motif.hpp"

MotifL2::MotifL2(bool transformed): MotifSobol(transformed) {};

std::variant<MotifPure::indexField,KMA::Mfield>
MotifL2::compute_motif(const arma::urowvec& v_dom,
                       const KMA::ivector& s_k,
                       const KMA::vector& p_k,
                       const KMA::Mfield& Y,
                       double m) const
{
  return compute_motif_helper<false>(v_dom,s_k,p_k,Y,m);
}
    
void MotifL2::set_parameters(const Parameters & newParameters){
    _transformed = newParameters._transformed;
}
    
void MotifL2::elongate_motifs(KMA::Mfield& V_new,
                              std::vector<arma::urowvec>& V_dom,
                              KMA::imatrix& S_k,const KMA::matrix& P_k,
                              const KMA::Mfield& Y,const KMA::matrix& D,
                              const Parameters& param,
                              const std::shared_ptr<PerformanceIndexAB>& perf,
                              const std::shared_ptr<Dissimilarity>& diss,
                              const Rcpp::Function & quantile_func) const
{
  return elongate_motifs_helper<false>(V_new,V_dom,S_k,P_k,Y,D,param,perf,diss,quantile_func);
}


