//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/util/f77.h>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>

using namespace std;

DensityFit::DensityFit(const int nbas, const int naux,
                                       const vector<shared_ptr<Atom> >& atoms,  const vector<vector<int> >& offsets,
                                       const vector<shared_ptr<Atom> >& aux_atoms,  const vector<vector<int> >& aux_offsets, const double throverlap)
   : nbasis_(nbas), naux_(naux) {

  // this will be distributed in the future.
  data_ = new double[nbasis_*nbasis_*naux];
  fill(data_, data_+nbasis_*nbasis_*naux, 0.0);

  vector<shared_ptr<Shell> > basis; 
  vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = offsets[cnt]; 
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }
  const int size = basis.size();

  // some info for auxiliary (i.e., DF) basis set
  vector<shared_ptr<Shell> > aux_basis; 
  vector<int> aux_offset;
  cnt = 0;
  for (auto aiter = aux_atoms.begin(); aiter != aux_atoms.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    aux_basis.insert(aux_basis.end(), tmp.begin(), tmp.end());  
    const vector<int> tmpoff = aux_offsets[cnt]; 
    aux_offset.insert(aux_offset.end(), tmpoff.begin(), tmpoff.end());
  }
  const int aux_size = aux_basis.size();

  const std::shared_ptr<Shell> b3(new Shell(basis.front()->spherical()));

  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<Shell>  b0 = basis[i0];
    const int b0offset = offset[i0]; 
    const int b0size = b0->nbasis();
    for (int i1 = i0; i1 != size; ++i1) {
      const shared_ptr<Shell>  b1 = basis[i1];
      const int b1offset = offset[i1]; 
      const int b1size = b1->nbasis();
      for (int i2 = 0; i2 != aux_size; ++i2) {
        const shared_ptr<Shell>  b2 = aux_basis[i2];
        const int b2offset = aux_offset[i2]; 
        const int b2size = b2->nbasis();

        vector<shared_ptr<Shell> > input;
        input.push_back(b3);
        input.push_back(b2);
        input.push_back(b1);
        input.push_back(b0);

        // TODO if I turn on primitive screening, it is broken.
        ERIBatch eribatch(input, 0.0);
        eribatch.compute();
        const double* eridata = eribatch.data();

        // all slot in
        for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
          for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {  
            for (int j2 = b2offset; j2 != b2offset + b2size; ++j2, ++eridata) {  
              data_[j2+naux*(j1+nbasis_*j0)] = data_[j2+naux*(j0+nbasis_*j1)] = *eridata;
            }
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  data2_ = new double[naux*naux];
  fill(data2_, data2_+naux*naux, 0.0);

  for (int i0 = 0; i0 != aux_size; ++i0) {
    const std::shared_ptr<Shell>  b0 = aux_basis[i0];
    const int b0offset = aux_offset[i0]; 
    const int b0size = b0->nbasis();
    for (int i1 = i0; i1 != aux_size; ++i1) {
      const std::shared_ptr<Shell>  b1 = aux_basis[i1];
      const int b1offset = aux_offset[i1]; 
      const int b1size = b1->nbasis();

      std::vector<std::shared_ptr<Shell> > input;
      input.push_back(b1);
      input.push_back(b3);
      input.push_back(b0);
      input.push_back(b3);

      // TODO if I turn on primitive screening, it is broken.
      ERIBatch eribatch(input, 0.0);
      eribatch.compute();
      const double* eridata = eribatch.data();

      for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {  
        for (int j1 = b1offset; j1 != b1offset + b1size; ++j1, ++eridata) {  
          data2_[j1+j0*naux] = data2_[j0+j1*naux] = *eridata; 
        }
      }
    }
  }

  const int lwork = 5*naux;
  double* vec = new double[naux];
  double* work = new double[std::max(lwork,naux*naux)];

  int info;
  dsyev_("V", "U", &naux, data2_, &naux, vec, work, &lwork, &info); 
  if (info) throw std::runtime_error("dsyev failed in DF fock builder");

  for (int i = 0; i != naux; ++i)
    vec[i] = vec[i] > throverlap ? 1.0/std::sqrt(std::sqrt(vec[i])) : 0.0;
  for (int i = 0; i != naux; ++i) {
    for (int j = 0; j != naux; ++j) {
      data2_[j+i*naux] *= vec[i];
    }
  }
  // work now contains -1/2
  dgemm_("N", "T", naux, naux, naux, 1.0, data2_, naux, data2_, naux, 0.0, work, naux); 
  copy(work, work+naux*naux, data2_);

  delete[] vec;
  delete[] work;
}


DensityFit::~DensityFit() {
  delete[] data_;
  delete[] data2_;
}

