//
// author : Toru Shiozaki
//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/scf/f77.h>
#include <src/fci/mofile.h>
#include <src/rysint/eribatch.h>
#include <src/scf/scf.h>

using namespace boost;
using namespace std;

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Shell> RefShell;

MOFile::MOFile(const shared_ptr<Geometry> geom, const shared_ptr<Coeff> cmo) : geom_(geom), coeff_(cmo) {
  {
    Filename tmpf;
    filename_ = tmpf.filename_next();
  }
  boost::shared_ptr<std::fstream> tmp(new std::fstream(filename_.c_str(), std::ios::out | std::ios::trunc | std::ios::binary));
  file_ = tmp;

  { // prepare offset and basis
    typedef boost::shared_ptr<Atom> RefAtom;
    typedef boost::shared_ptr<Shell> RefShell;

    const std::vector<RefAtom> atoms = geom_->atoms();
    int cnt = 0;
    for (std::vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
      const std::vector<RefShell> tmp = (*aiter)->shells();
      basis_.insert(basis_.end(), tmp.begin(), tmp.end());  
      const std::vector<int> tmpoff = geom_->offset(cnt); 
      offset_.insert(offset_.end(), tmpoff.begin(), tmpoff.end());
    }
  }
}

MOFile::~MOFile() {
  unlink(filename_.c_str());
}

// I don't care the efficiency at all!!
// >>>>>>>>>>>>  ALL INCORE!!! <<<<<<<<<<<
void MOFile::create_Jiiii(const int nstart, const int nfence) {
  // first compute all the AO integrals in core

  const int nocc = nfence - nstart;
  const int size = basis_.size(); // number of shells
  const int nbasis = geom_->nbasis();
  const size_t aointsize = nbasis*nbasis*nbasis*nbasis; 
  double* aobuff = new double[aointsize];
  double* first = new double[nbasis*nbasis*nbasis*nocc];

  // some stuffs for blas
  const int unit = 1;
  const double one = 1.0;
  const double zero = 0.0;
  double* cdata = coeff_->data() + nstart*nbasis;

  cout << "  - AO integrals are computed and stored in core" << endl;

  // one electron part
  Hcore hcore(geom_); 
  dgemm_("n","n",&nbasis,&nocc,&nbasis,&one,hcore.data(),&nbasis,cdata,&nbasis,&zero,aobuff,&nbasis);
  mo1e_.resize(nocc*nocc);
  dgemm_("t","n",&nocc,&nocc,&nbasis,&one,cdata,&nbasis,aobuff,&nbasis,&zero,&mo1e_[0],&nocc);

  for (int i0 = 0; i0 != size; ++i0) {
    const int b0offset = offset_[i0]; 
    const int b0size = basis_[i0]->nbasis();
    for (int i1 = 0; i1 != size; ++i1) {
      const int b1offset = offset_[i1];
      const int b1size = basis_[i1]->nbasis();
      for (int i2 = 0; i2 != size; ++i2) {
        const int b2offset = offset_[i2];
        const int b2size = basis_[i2]->nbasis();
        for (int i3 = 0; i3 != size; ++i3) {
          const int b3offset = offset_[i3]; 
          const int b3size = basis_[i3]->nbasis();
          vector<RefShell> input;
          input.push_back(basis_[i3]);
          input.push_back(basis_[i2]);
          input.push_back(basis_[i1]);
          input.push_back(basis_[i0]);

          ERIBatch eribatch(input, 1.0);
          eribatch.compute();
          
          const double* eridata = eribatch.data();
          // what a bad code!
          int ioff = 0;
          for (int i = b0offset; i != b0offset+b0size; ++i) {
            for (int j = b1offset; j != b1offset+b1size; ++j) {
              for (int k = b2offset; k != b2offset+b2size; ++k) {
                for (int l = b3offset; l != b3offset+b3size; ++l, ++ioff) {
                  aobuff[l+nbasis*(k+nbasis*(j+nbasis*i))] = eridata[ioff];
                }
              }
            }
          } 
          
        }
      }
    }
  }

  const int n = aointsize;
//cout << ddot_(&n, aobuff, &unit, aobuff, &unit) << endl; 
  const int nn = nbasis*nbasis;
  const int nnn = nbasis*nbasis*nbasis;
  dgemm_("n","n",&nnn,&nocc,&nbasis, &one, aobuff,&nnn,cdata,&nbasis, &zero,first,&nnn);

  for (int i = 0; i != nocc; ++i) {
    dgemm_("n","n",&nn,&nocc,&nbasis, &one, first+nnn*i,&nn,cdata,&nbasis, &zero,aobuff+nn*nocc*i,&nn);
  }

  const int mm = nocc*nocc;
  mytranspose_(aobuff,&nn,&mm,first);

  const int nmm = nbasis * mm; 
  dgemm_("n","n",&nmm,&nocc,&nbasis,&one,first,&nmm,cdata,&nbasis,&zero,aobuff,&nmm);

  for (int i = 0; i != nocc; ++i) {
    dgemm_("n","n",&mm,&nocc,&nbasis,&one,aobuff+nmm*i,&mm,cdata,&nbasis,&zero,first+mm*nocc*i,&mm);
  }
  delete[] aobuff;
  mo2e_.resize(mm*mm);
  copy(first,first+mm*mm,&mo2e_[0]);
  delete[] first;
}
