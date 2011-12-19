//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_hcore_h
#define __src_scf_hcore_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <memory>

class Hcore : public Matrix1e {
  protected:
    void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int);

  public:
    Hcore(const std::shared_ptr<Geometry>);
    ~Hcore();

};

#endif

