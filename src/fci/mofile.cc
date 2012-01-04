//
// author : Toru Shiozaki
//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <src/util/f77.h>
#include <src/fci/mofile.h>
#include <src/rysint/eribatch.h>
#include <src/scf/scf.h>

using namespace std;

typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

MOFile::MOFile(const shared_ptr<Geometry> geom, const shared_ptr<Reference> ref) : geom_(geom), ref_(ref) {

  do_df_ = geom->df().get();
do_df_ = false;

  typedef std::shared_ptr<Atom> RefAtom;
  typedef std::shared_ptr<Shell> RefShell;

  const std::vector<RefAtom> atoms = geom_->atoms();
  int cnt = 0;
  for (std::vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const std::vector<RefShell> tmp = (*aiter)->shells();
    basis_.insert(basis_.end(), tmp.begin(), tmp.end());  
    const std::vector<int> tmpoff = geom_->offset(cnt); 
    offset_.insert(offset_.end(), tmpoff.begin(), tmpoff.end());
  }
}

MOFile::~MOFile() {
}



double MOFile::create_Jiiii(const int nstart, const int nfence) {

  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;
  nbasis_ = geom_->nbasis();
  const int nbasis = nbasis_;
  const int nocc = nocc_;
  const int size = basis_.size(); // number of shells
  const size_t aointsize = nbasis*nbasis*nbasis*nbasis; 
  double *aobuff, *first;
  if (!do_df_) {
    aobuff = new double[aointsize];
    first = new double[nbasis*nbasis*nbasis*nocc];
    cout << "  - AO integrals are computed and stored in core" << endl << endl;
  } else {
    first = new double[nocc*nocc*nocc*nocc];
    aobuff = new double[nbasis*nocc];
  }

  // some stuffs for blas
  const int unit = 1;
  const double one = 1.0;
  const double zero = 0.0;
  double* cdata = ref_->coeff()->data() + nstart*nbasis;
  const int mm = nocc*nocc;

  // one electron part
  double core_energy = 0.0;
  {
    const bool df_ = do_df_;
    if (df_) {
      shared_ptr<Fock<1> > fock0(new Fock<1>(geom_, ref_->hcore()));
      if (nstart != 0) {
        shared_ptr<Matrix1e> den(new Matrix1e(ref_->coeff()->form_core_density_rhf()));
        shared_ptr<Fock<1> > fock1(new Fock<1>(geom_, fock0, den, ref_->schwarz()));
        core_energy = (*den * (*ref_->hcore()+*fock1)).trace();
        fock0 = fock1;
      }
      fock0->symmetrize();
      dgemm_("n","n",&nbasis,&nocc,&nbasis,&one,fock0->data(),&nbasis,cdata,&nbasis,&zero,aobuff,&nbasis);
    } else {
      shared_ptr<Fock<0> > fock0(new Fock<0>(geom_, ref_->hcore()));
      if (nstart != 0) {
        shared_ptr<Matrix1e> den(new Matrix1e(ref_->coeff()->form_core_density_rhf()));
        shared_ptr<Fock<0> > fock1(new Fock<0>(geom_, fock0, den, ref_->schwarz()));
        core_energy = (*den * (*ref_->hcore()+*fock1)).trace();
        fock0 = fock1;
      }
      fock0->symmetrize();
      dgemm_("n","n",&nbasis,&nocc,&nbasis,&one,fock0->data(),&nbasis,cdata,&nbasis,&zero,aobuff,&nbasis);
    }
  }
  mo1e_.resize(nocc*nocc);
  dgemm_("t","n",&nocc,&nocc,&nbasis,&one,cdata,&nbasis,aobuff,&nbasis,&zero,&mo1e_[0],&nocc);

  /////////////// non Df builder ///////////////
  if (!do_df_) {
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
            for (int i = b0offset; i != b0offset+b0size; ++i)
              for (int j = b1offset; j != b1offset+b1size; ++j)
                for (int k = b2offset; k != b2offset+b2size; ++k)
                  for (int l = b3offset; l != b3offset+b3size; ++l, ++ioff)
                    aobuff[l+nbasis*(k+nbasis*(j+nbasis*i))] = eridata[ioff];
          }
        }
      }
    }

    const int n = aointsize;
    const int nn = nbasis*nbasis;
    const int nnn = nbasis*nbasis*nbasis;
    dgemm_("n","n",&nnn,&nocc,&nbasis, &one, aobuff,&nnn,cdata,&nbasis, &zero,first,&nnn);

    for (int i = 0; i != nocc; ++i)
      dgemm_("n","n",&nn,&nocc,&nbasis, &one, first+nnn*i,&nn,cdata,&nbasis, &zero,aobuff+nn*nocc*i,&nn);

    mytranspose_(aobuff,&nn,&mm,first);

    for (int i = 0; i != nbasis; ++i)
      dgemm_("n","n",&mm,&nocc,&nbasis,&one,first+mm*nbasis*i,&mm,cdata,&nbasis,&zero,aobuff+mm*nocc*i,&mm);

    // aobuff here contains (rs|tx) with r running fastest. x: AO
    // Anyway this is not the way. DF will be used. TODO
    mo2e_1ext_.resize(mm*nbasis*nocc);
    copy(aobuff, aobuff+mm*nbasis*nocc, mo2e_1ext_ptr());

    const int nmm = nocc * mm;
    dgemm_("n","n",&nmm,&nocc,&nbasis,&one,aobuff,&nmm,cdata,&nbasis,&zero,first,&nmm);

    delete[] aobuff;

  //////////////////// df builder ////////////////////////
  } else {
    delete[] aobuff;

    shared_ptr<DensityFit> dff = geom_->df();
    const int naux = dff->naux();
    // we want to store half-transformed quantity for latter convenience
    mo2e_1ext_.resize(nocc*naux*nbasis);

    // first half transformation
    for (size_t i = 0; i != nbasis_; ++i)
      dgemm_("N", "N", naux, nocc, nbasis, 1.0, dff->data_3index()+i*nbasis*naux, naux, cdata, nbasis, 0.0, &(mo2e_1ext_[i*naux*nocc]), naux);

    double* buf = new double[naux*nocc*nocc];
    double* buf2 = new double[naux*nocc*nocc];
    dgemm_("n", "n", naux*nocc, nocc, nbasis, 1.0, &(mo2e_1ext_[0]), naux*nocc, cdata, nbasis, 0.0, buf2, naux*nocc);
    dgemm_("n", "n", naux, nocc*nocc, naux, 1.0, dff->data_2index(), naux, buf2, naux, 0.0, buf, naux);
    dgemm_("t", "n", nocc*nocc, nocc*nocc, naux, 1.0, buf, naux, buf, naux, 0.0, first, nocc*nocc);
    delete[] buf;
    delete[] buf2;
  }


  // mo2e is compressed
  sizeij_ = nocc*(nocc+1)/2;
  mo2e_.resize(sizeij_*sizeij_);
  // TODO very bad
  int ijkl = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j <= i; ++j) {
      const int ijo = (j + i*nocc)*nocc*nocc;
      for (int k = 0; k != nocc; ++k) {
        for (int l = 0; l <= k; ++l, ++ijkl) {
          mo2e_[ijkl] = first[l+k*nocc+ijo]; 
        }
      }
    }
  }

  // h'kl = hkl - 0.5 sum_j (kj|jl)
  vector<double> buf3(sizeij_);
  int ij = 0;
  for (int i=0; i!=nocc; ++i) {
    for (int j=0; j<=i; ++j, ++ij) {
      buf3[ij] = mo1e_[j+i*nocc];
      for (int k=0; k!=nocc; ++k) {
        buf3[ij] -= 0.5*first[(k+j*nocc)*mm+(k+i*nocc)];
      }
    }
  }
  copy(buf3.begin(), buf3.end(), mo1e_.begin());
  delete[] first;
  return core_energy; // TODO this is not a good way of implementation...
}


void MOFile::update_1ext_ints(const vector<double>& coeff) {
  // in the case of no DF
  double* buf = new double[mo2e_1ext_.size()];
  if (!do_df_) {
    const double* start = &(coeff[0]);
    double* data_ = &(mo2e_1ext_[0]);
    const int n2 = nocc_*nocc_;
    const int nb = nocc_*nbasis_;
    assert(n2*nb == mo2e_1ext_.size());

    mytranspose_(data_, &n2, &nb, buf);
    // now buf=(nx|nn).
    dgemm_("N", "N", nb*nocc_, nocc_, nocc_, 1.0, buf, nb*nocc_, start, nocc_, 0.0, data_, nb*nocc_); 
    for (int i = 0; i != nocc_; ++i)
      dgemm_("N", "N", nb, nocc_, nocc_, 1.0, data_+i*nb*nocc_, nb, start, nocc_, 0.0, buf+i*nb*nocc_, nb); 
    dgemm_("T", "N", nocc_, nb*nocc_, nocc_, 1.0, start, nocc_, buf, nocc_, 0.0, data_, nocc_); 
    // slot in place
    mytranspose_(data_, &nb, &n2, buf);
  } else {
    // half transformed DF is rotated.
    const int naux = geom_->df()->naux();
    for (int i = 0; i != nbasis_; ++i) 
      dgemm_("N", "N", naux, nocc_, nocc_, 1.0, &(mo2e_1ext_[i*naux*nocc_]), naux, &(coeff[0]), nocc_, 0.0, buf+i*naux*nocc_, naux); 
  }
  copy(buf, buf+mo2e_1ext_.size(), &(mo2e_1ext_[0]));
  delete[] buf; 
}
