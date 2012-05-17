//
// Author : Toru Shiozaki
// Date   : May 2012
//

#ifndef __SRC_PROP_DIPOLE_H
#define __SRC_PROP_DIPOLE_H

#include <memory>
#include <src/scf/matrix1e.h>

class Dipole {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Matrix1e> den_;

  public:
    Dipole(std::shared_ptr<const Geometry>, std::shared_ptr<const Matrix1e>);
    ~Dipole();

    void compute();
};

#endif
