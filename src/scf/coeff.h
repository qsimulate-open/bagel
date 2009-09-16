//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_coeff_h
#define __src_scf_coeff_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>

class Coeff : public Matrix1e {
  protected:

  public:
    Coeff(const Matrix1e&);
    ~Coeff();

    Matrix1e form_density_rhf() const;
};

#endif
