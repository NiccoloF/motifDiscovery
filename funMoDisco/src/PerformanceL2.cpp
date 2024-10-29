#include "PerformanceIndex.hpp"

double PerformanceL2::compute_Jk(const KMA::Mfield& V,const KMA::ivector& s_k,
                                 const KMA::vector& p_k,const KMA::Mfield& Y,
                                 const KMA::vector& w,int m,double c_k, 
                                 const KMA::uvector& keep_k,
                                 const std::shared_ptr<Dissimilarity>& diss) const
{
  return compute_Jk_helper<false>(V,s_k,p_k,Y,w,m,c_k,keep_k,diss);
}
