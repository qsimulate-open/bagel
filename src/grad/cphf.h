//
// Author : Toru Shiozaki
// Date   : May 2012
//

#ifndef __SRC_GRAD_CPHF_H
#define __SRC_GRAD_CPHF_H

#include <src/scf/matrix1e.h>
#include <memory>
#include <src/util/linear.h>
#include <src/wfn/reference.h>

class CPHF {
  protected:
    const std::shared_ptr<Linear<Matrix1e> > solver_;
    const std::shared_ptr<const Matrix1e> grad_;
    const std::vector<double> eig_;
    const std::shared_ptr<DF_Half> halfjj_;
    const std::shared_ptr<const Reference> ref_;
    const std::shared_ptr<const Geometry> geom_;
    const int ncore_;

  public:
    CPHF(const std::shared_ptr<const Matrix1e> grad, const std::vector<double>& eig,
         const std::shared_ptr<DF_Half> half, const std::shared_ptr<const Reference> g, const int ncore);
    ~CPHF() {};

    std::shared_ptr<Matrix1e> solve() const;

};

#endif

