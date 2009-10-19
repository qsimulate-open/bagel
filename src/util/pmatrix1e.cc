//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <src/util/pmatrix1e.h>
#include <src/pscf/f77.h>
#include <src/macros.h>

typedef std::complex<double> Complex;

using namespace std;


// I have supposed that ndim_ is the leading dimension in the data.

PMatrix1e::PMatrix1e(const boost::shared_ptr<PGeometry> pg) 
: geom_(pg), nbasis_(pg->nbasis()),
  blocksize_(pg->nbasis() * pg->nbasis()), totalsize_((2 * pg->K() + 1) * pg->nbasis() * pg->nbasis()) {

  mdim_ = ndim_ = nbasis_;
  boost::shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp; 
}


PMatrix1e::PMatrix1e(const boost::shared_ptr<PGeometry> pg, const int ldn, const int ldm)
: geom_(pg), nbasis_(pg->nbasis()),
  blocksize_(ldn * ldm), totalsize_((2 * pg->K() + 1) * ldn * ldm) {

  mdim_ = ldm;
  ndim_ = ldn;
  boost::shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
}


// Changes the leading dimension, filling zero.
PMatrix1e::PMatrix1e(const boost::shared_ptr<PMatrix1e> source, const int ldn, const int ldm)
: geom_(source->geom()), nbasis_(source->nbasis()),
  blocksize_(ldn * ldm), totalsize_((2 * source->K() + 1) * ldn * ldm) {

  assert(ldn >= nbasis_);

  ndim_ = ldn;
  mdim_ = ldm;
  boost::shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
  const boost::shared_ptr<PData> sdata = source->data();
  const int ld_source = source->ndim();
  const int sblocksize = source->blocksize();

  for (int i = -K(); i != max(K(), 1); ++i) {
    complex<double>* current = data_->front() + (i + K()) * blocksize_;
    complex<double>* csource = sdata->front() + (i + K()) * sblocksize;

    for (int j = 0; j != mdim_; ++j, current+=ldn, csource+=ld_source) {
      copy(csource, csource+ld_source, current);
    }
  }
}


// Changes number of columns; removes first ncut columns.
PMatrix1e::PMatrix1e(const boost::shared_ptr<PMatrix1e> source, const int mcut)
: geom_(source->geom()), nbasis_(source->nbasis()),
  blocksize_(source->ndim() * (source->mdim()-mcut)),
  totalsize_((2 * source->K() + 1) * source->ndim() * (source->mdim()-mcut)) {

  assert(mcut < source->mdim());

  ndim_ = source->ndim();
  mdim_ = source->mdim()-mcut;
  boost::shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
  const boost::shared_ptr<PData> sdata = source->data();
  const int sblocksize = source->blocksize();

  for (int i = -K(); i != max(K(), 1); ++i) {
    complex<double>* current = data_->front() + (i + K()) * blocksize_;
    complex<double>* csource = sdata->front() + (i + K()) * sblocksize;

    copy(csource+ndim_*mcut, csource+ndim_*(mdim_+mcut), current);
  }
}



PMatrix1e::~PMatrix1e() {
}


PMatrix1e PMatrix1e::operator+(const PMatrix1e& add) const {
  assert(add.totalsize() == totalsize_);
  PMatrix1e out(geom_);
  for (int j = 0; j != totalsize_; ++j) (*out.data_)[j] = (*add.data_)[j] + (*data_)[j]; 
  return out;
}


PMatrix1e PMatrix1e::operator-(const PMatrix1e& sub) const {
  assert(sub.totalsize() == totalsize_);
  PMatrix1e out(geom_);
  #pragma omp parallel for
  for (int j = 0; j < totalsize_; ++j) (*out.data_)[j] = (*data_)[j] - (*sub.data_)[j]; 
  return out;
}


PMatrix1e& PMatrix1e::operator=(const PMatrix1e& source) {
  assert(source.totalsize() == totalsize_);
  ::memcpy(data_->front(), source.data_->front(), totalsize_ * sizeof(Complex)); 
  return *this; 
}


PMatrix1e& PMatrix1e::operator+=(const PMatrix1e& source) {
  assert(source.totalsize() == totalsize_);
  const int unit = 1;
  const Complex one(1.0, 0.0);
  zaxpy_(&totalsize_, &one, source.data_->front(), &unit, data_->front(), &unit);
  return *this; 
}


PMatrix1e PMatrix1e::ft() const {
  PMatrix1e out(geom_);

  const int unit = 1;
  for (int k = -K(); k < max(K(), 1); ++k) {
    const Complex kq(0.0, K() != 0 ? (DPI * k) / K() : 0.0);
    const int boffset = (k + K()) * blocksize_;
    for (int m = -S(), mcount = 0; m <= S(); ++m, ++mcount) {
      const int moffset = mcount * blocksize_;
      Complex factor = exp(- kq * static_cast<double>(m)); // caution 
      zaxpy_(&blocksize_, &factor, data_->front() + moffset, &unit, out.data_->front() + boffset, &unit);
    }
  }
  return out;
}


PMatrix1e PMatrix1e::bft() const {
  PMatrix1e out(geom_);

  const int unit = 1;
  if (K() != 0) {
    for (int m = -K(); m <= K(); ++m) {
      const int mcount = m + K();
      const int moffset = mcount * blocksize_;
      for (int k = -K(), kcount = 0; k != K(); ++k, ++kcount) {
        const Complex kq(0.0, (DPI * k) / K());
        const int koffset = kcount * blocksize_;
        Complex factor = exp(kq * static_cast<double>(m)) / (2.0 * K());
        zaxpy_(&blocksize_, &factor, data_->front() + koffset, &unit, out.data_->front() + moffset, &unit);
      }
    }
  } else {
    zcopy_(&totalsize_, data_->front(), &unit, out.data_->front(), &unit);
  }
  return out;
}


void PMatrix1e::init() {
  typedef boost::shared_ptr<Atom> RefAtom;
  typedef boost::shared_ptr<Shell> RefShell;

  const std::vector<RefAtom> atoms = geom_->atoms();

  const std::vector<std::vector<int> > offsets = geom_->offsets();
  const int natom = geom_->natom();
  //#pragma omp parallel for
  for (int iatom0 = 0; iatom0 < natom; ++iatom0) {
    const RefAtom catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const std::vector<int> coffset0 = offsets[iatom0];
    const std::vector<RefShell> shell0 = catom0->shells();
    for (int m = -K(), mcount = 0; m <= K(); ++m, ++mcount) {
      const int blockoffset = mcount * blocksize_;
      for (int iatom1 = 0; iatom1 != geom_->natom(); ++iatom1) {
        double disp[3];
        disp[0] = disp[1] = 0.0; 
        disp[2] = geom_->A() * m;
        const RefAtom catom1(new Atom(*(atoms[iatom1]), disp));
        const int numshell1 = catom1->shells().size();
        const std::vector<int> coffset1 = offsets[iatom1];
        const std::vector<RefShell> shell1 = catom1->shells();
  
        for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
          const int offset0 = coffset0[ibatch0]; 
          RefShell b0 = shell0[ibatch0];
          for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
            const int offset1 = coffset1[ibatch1];
            RefShell b1 = shell1[ibatch1];
            std::vector<RefShell> input;
            input.push_back(b1);
            input.push_back(b0);
  
            computebatch(input, offset0, offset1, nbasis_, blockoffset);
  
          }
        } 
      }
    } 
  }
};


PMatrix1e PMatrix1e::operator*(const PMatrix1e& o) const {
  const int unit = 1;
  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();

  PMatrix1e out(geom_, l, n);
  const boost::shared_ptr<PData> odata = o.data();
  boost::shared_ptr<PData> outdata = out.data();

  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset1 = icount * blocksize_;
    const int boffset2 = icount * o.blocksize();
    const int boffset3 = icount * out.blocksize();
    const Complex* cdata = data_->pointer(boffset1);
    const Complex* codata = odata->pointer(boffset2);
    Complex* coutdata = outdata->pointer(boffset3);
    zgemm_("N", "N", &l, &n, &m, &one, cdata, &ndim_, codata, &mdim_, &zero, coutdata, &ndim_);
  }

  return out;
};


PMatrix1e PMatrix1e::operator%(const PMatrix1e& o) const {

  const int unit = 1;
  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);
  const int l = mdim_;
  const int m = ndim_; 
  assert(ndim_ == o.ndim());
  const int n = o.mdim(); 
  const boost::shared_ptr<PData> odata = o.data();

  PMatrix1e out(geom_, l, n);
  boost::shared_ptr<PData> outdata = out.data();

  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset1 = icount * blocksize_;
    const int boffset2 = icount * o.blocksize();
    const int boffset3 = icount * out.blocksize();
    const Complex* cdata = data_->pointer(boffset1);
    const Complex* codata = odata->pointer(boffset2);
    Complex* coutdata = outdata->pointer(boffset3);
    zgemm_("C", "N", &l, &n, &m, &one, cdata, &ndim_, codata, &ndim_, &zero, coutdata, &mdim_);
  }

  return out;
}


void PMatrix1e::svd(boost::shared_ptr<PMatrix1e> U, boost::shared_ptr<PMatrix1e> V) {
  assert(U->ndim() == ndim_ && U->mdim() == ndim_);
  assert(V->ndim() == mdim_ && V->mdim() == mdim_);
  const int lwork = 5 * min(ndim_, mdim_) + max(ndim_, mdim_);
  complex<double>* work = new complex<double>[lwork];
  double* rwork = new double[lwork];
  double* S = new double[min(ndim_, mdim_)];
/*
  SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, RWORK, INFO )
 */
  for (int i = -K(); i < max(K(), 1); ++i) {
    complex<double>* cblock = data_->pointer(blocksize_ * (i + K()));
    complex<double>* ublock = U->data()->pointer(U->blocksize() * (i + K()));
    complex<double>* vblock = V->data()->pointer(V->blocksize() * (i + K()));
    int info = 0;
    zgesvd_("A", "A", &ndim_, &mdim_, cblock, &ndim_, S, ublock, &ndim_, vblock, &mdim_, work, &lwork, rwork, &info);
    assert(info == 0);
  }
  delete[] S;
  delete[] work;
  delete[] rwork;
}


void PMatrix1e::diagonalize(double* eig) {
  assert(ndim_ == mdim_);
  const int lwork = nbasis_ * 6;

  #pragma omp parallel for
  for (int i = -K(); i <= K(); ++i) {
    Complex* work = new Complex[lwork];
    double* rwork = new double[nbasis_ * 3 - 2];
    const int icount = i + K();
    const int boffset = icount * blocksize_;
    const int eoffset = icount * nbasis_;
    Complex* cdata = data_->pointer(boffset);
    double* ceig = &eig[eoffset];

    int info_diagonalize = 0;
    zheev_("V", "L", &ndim_, cdata, &nbasis_, ceig, work, &lwork, rwork, &info_diagonalize); 
    assert(info_diagonalize == 0);
    delete[] work;
    delete[] rwork;
  }

}


void PMatrix1e::print() const {
  int index = 0;
  for (int k = -K(); k <= K(); ++k) {
    cout << "K = " << setw(2) << k << ":" << endl;
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j, ++index) {
        cout << setw(14) << setprecision(6) << (*data_)[index] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}


void PMatrix1e::rprint() const {
  int index = 0;
  for (int k = -K(); k <= K(); ++k) {
    cout << "K = " << setw(2) << k << ":" << endl;
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j, ++index) {
        cout << setw(10) << setprecision(6) << (*data_)[index].real() * 2.0 << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}


void PMatrix1e::hermite() {
  #pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    const int kcount = k + K();
    const int koffset = kcount * blocksize_; 
    Complex* dat = data_->pointer(koffset);
    for (int i = 0; i != nbasis_; ++i) {
      for (int j = 0; j != nbasis_; ++j) {
        const Complex a = dat[j + i * nbasis_]; 
        const Complex b = dat[i + j * nbasis_]; 
        dat[j + i * nbasis_] = (a + conj(b)) * 0.5;
        dat[i + j * nbasis_] = (b + conj(a)) * 0.5;
      }
    }
  }
}


void PMatrix1e::real() {
  #pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    const int kcount = k + K();
    const int koffset = kcount * blocksize_; 
    Complex* dat = data_->pointer(koffset);
    for (int i = 0; i != nbasis_; ++i) {
      for (int j = 0; j != nbasis_; ++j) {
        dat[j + i * nbasis_] = (dat[j + i * nbasis_]).real();
      }
    }
  }
}


void PMatrix1e::scale(const Complex a) {
  const int unit = 1;
  zscal_(&totalsize_, &a, data_->front(), &unit);
}


void PMatrix1e::zaxpy(const Complex a, const PMatrix1e& o) {
  const int unit = 1;
  const Complex* odata = o.data()->front();
  zaxpy_(&totalsize_, &a, odata, &unit, data_->front(), &unit); 
}


void PMatrix1e::zaxpy(const Complex a, const boost::shared_ptr<PMatrix1e> o) {
  const int unit = 1;
  const Complex* odata = o->data()->front();
  zaxpy_(&totalsize_, &a, odata, &unit, data_->front(), &unit); 
}


const Complex PMatrix1e::zdotc(const PMatrix1e& o) const {
  assert(o.totalsize() == totalsize_);
  const int unit = 1;
  const Complex* odata = o.data()->front();
  Complex out;
#ifndef ZDOT_RETURN
  zdotc_(&out, &totalsize_, data_->front(), &unit, odata, &unit); 
#else
  out = zdotc_(&totalsize_, data_->front(), &unit, odata, &unit); 
#endif
  return out;
}


const Complex PMatrix1e::zdotc(const boost::shared_ptr<PMatrix1e> o) const {
  assert(o->totalsize() == totalsize_);
  const int unit = 1;
  const Complex* odata = o->data()->front();
  Complex out;
#ifndef ZDOT_RETURN
  zdotc_(&out, &totalsize_, data_->front(), &unit, odata, &unit); 
#else
  out = zdotc_(&totalsize_, data_->front(), &unit, odata, &unit); 
#endif
  return out;
}


const double PMatrix1e::rms() const {
  const int unit = 1;
  Complex zdot;
#ifndef ZDOT_RETURN
  zdotc_(&zdot, &totalsize_, data_->front(), &unit, data_->front(), &unit); 
#else
  zdot = zdotc_(&totalsize_, data_->front(), &unit, data_->front(), &unit); 
#endif
  return ::sqrt(abs(zdot) / totalsize_);
}


const double PMatrix1e::trace() const {
  double out = 0.0;
  int scount = 0;
  for (int s = -L(); s <= L(); ++s, ++scount) {
    const int boffset =  scount * blocksize_; 
    for (int i = 0; i != nbasis_; ++i)
      out += ((*data_)[boffset + i * nbasis_ + i]).real(); 
  }
  return out;
}


