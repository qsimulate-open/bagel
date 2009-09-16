//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pcoeff_h
#define __src_pscf_pcoeff_h

#include <src/util/pmatrix1e.h>

class PCoeff : public PMatrix1e {
  protected:

  public:
    PCoeff(const PMatrix1e&);
    ~PCoeff();

    PMatrix1e form_density_rhf();
};

#endif

