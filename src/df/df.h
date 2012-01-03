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
    double* data_;
    double* data2_;
    size_t nbasis_;
    size_t naux_;

  public:
    DensityFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr);
    ~DensityFit();

    const double* const data_3index() const { return data_; };
    const double* const data_2index() const { return data2_; };

    size_t nbasis() const { return nbasis_; };
    size_t naux() const { return naux_; };

}; 

#endif

