//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_fock_base_h
#define __src_scf_fock_base_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <src/scf/hcore.h>
#include <memory>
#include <vector>

class Fock_base : public Matrix1e {
  protected:
    const std::shared_ptr<Fock_base> previous_;
    const std::shared_ptr<Matrix1e> density_;
    void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int);

    // virtual function that is to be defined in the derived class
    virtual void fock_two_electron_part() = 0;
    void fock_one_electron_part();

//  for debug
//  void slater_two_electron_part();

    std::vector<double> schwarz_;
    double schwarz_thresh_;

  public:
    Fock_base(const std::shared_ptr<Geometry>, const std::shared_ptr<Fock_base>, const std::shared_ptr<Matrix1e>, const std::vector<double>&);
    Fock_base(const std::shared_ptr<Geometry>, const std::shared_ptr<Hcore>);

    ~Fock_base();

};

#endif
