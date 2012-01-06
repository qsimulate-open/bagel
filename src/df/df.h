//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_DF_DensityFit_H
#define __NEWINT_DF_DensityFit_H

#include <vector>
#include <memory>
#include <src/scf/atom.h>

class DensityFit {
  protected:
    std::unique_ptr<double[]> data_;
    std::unique_ptr<double[]> data2_;
    size_t nbasis_;
    size_t naux_;

    double* data() { return data_.get(); };
    double* data2() { return data2_.get(); };
    const double* const data() const { return data_.get(); };
    const double* const data2() const { return data2_.get(); };

  public:
    DensityFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr);
    ~DensityFit() {};

    const double* const data_3index() const { return data(); };
    const double* const data_2index() const { return data2(); };

    size_t nbasis() const { return nbasis_; };
    size_t naux() const { return naux_; };

}; 

#endif

