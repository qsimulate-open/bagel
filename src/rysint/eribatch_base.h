//
// Author : Toru Shiozaki
// Date   : May 2012
//

#ifndef __SRC_RYSINT_ERIBATCH_BASE_H
#define __SRC_RYSINT_ERIBATCH_BASE_H

#include <src/rysint/macros.h>
#include <vector>
#include <src/scf/shell.h>
#include <src/rysint/rysint.h>
#include <tuple>

class ERIBatch_base : public RysInt{
  protected:
    void root_weight(const int ps); 
    void compute_ssss(const double);

  public:
    ERIBatch_base(const std::vector<std::shared_ptr<Shell> >& o, const double max_density, const int deriv) : RysInt(o) {
      const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
      deriv_rank_ = deriv;

      // determins if we want to swap shells
      set_swap_info(true);

      // stores AB and CD
      set_ab_cd();

      // set primsize_ and contsize_, as well as relevant members
      set_prim_contsizes();

      // sets angular info
      int asize_final, csize_final, asize_final_sph, csize_final_sph;
      std::tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = set_angular_info();

      // allocate
      allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);
      allocate_arrays(primsize_);

      compute_ssss(integral_thresh);

      root_weight(primsize_);

    };
    ~ERIBatch_base() {};

};

#endif
