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
    PCoeff(const std::shared_ptr<PGeometry>, const int ldn, const int ldm);
    ~PCoeff();

    PMatrix1e form_density_rhf(const bool return_ao=true) const;

    std::pair<std::shared_ptr<PCoeff>, std::shared_ptr<PCoeff> >
          split(const int nrow1, const int nrow2);

};

#endif

