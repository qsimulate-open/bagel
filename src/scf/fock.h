//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_fock_h
#define __src_scf_fock_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <src/scf/hcore.h>
#include <boost/shared_ptr.hpp>
#include <vector>

class Fock : public Matrix1e {
  protected:
    const boost::shared_ptr<Fock> previous_;
    const boost::shared_ptr<Matrix1e> density_;
    void computebatch(const std::vector<boost::shared_ptr<Shell> >&, const int, const int, const int);

    void fock_two_electron_part();
// for debug
    void slater_two_electron_part();
    std::vector<double> shwarz_;

  public:
    Fock(const boost::shared_ptr<Geometry>, const boost::shared_ptr<Fock>, const boost::shared_ptr<Matrix1e>, const std::vector<double>&);
    Fock(const boost::shared_ptr<Geometry>, const boost::shared_ptr<Hcore>);

    ~Fock();

};

#endif
