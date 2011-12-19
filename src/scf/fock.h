//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_fock_h
#define __src_scf_fock_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <src/scf/hcore.h>
#include <memory>
#include <vector>

class Fock : public Matrix1e {
  protected:
    const std::shared_ptr<Fock> previous_;
    const std::shared_ptr<Matrix1e> density_;
    void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int);

    void fock_two_electron_part();
// for debug
    void slater_two_electron_part();
    std::vector<double> shwarz_;

  public:
    Fock(const std::shared_ptr<Geometry>, const std::shared_ptr<Fock>, const std::shared_ptr<Matrix1e>, const std::vector<double>&);
    Fock(const std::shared_ptr<Geometry>, const std::shared_ptr<Hcore>);

    ~Fock();

};

#endif
