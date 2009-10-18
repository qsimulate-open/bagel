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

PMatrix1e::PMatrix1e(const boost::shared_ptr<PGeometry> pg) 
: geom_(pg), nbasis_(pg->nbasis()),
  blocksize_(pg->nbasis() * pg->nbasis()), totalsize_((2 * pg->K() + 1) * pg->nbasis() * pg->nbasis()) {

  mdim_ = ndim_ = nbasis_;
  boost::shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp; 
}


PMatrix1e::~PMatrix1e() {
}


PMatrix1e PMatrix1e::operator+(const PMatrix1e& add) const {
  PMatrix1e out(geom_);
  for (int j = 0; j != totalsize_; ++j) (*out.data_)[j] = (*add.data_)[j] + (*data_)[j]; 
  return out;
}


PMatrix1e PMatrix1e::operator-(const PMatrix1e& sub) const {
  PMatrix1e out(geom_);
  #pragma omp parallel for
  for (int j = 0; j < totalsize_; ++j) (*out.data_)[j] = (*data_)[j] - (*sub.data_)[j]; 
  return out;
}


PMatrix1e& PMatrix1e::operator=(const PMatrix1e& source) {
  ::memcpy(data_->front(), source.data_->front(), totalsize_ * sizeof(Complex)); 
  return *this; 
}


PMatrix1e& PMatrix1e::operator+=(const PMatrix1e& source) {
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

  PMatrix1e out(geom_);
  const boost::shared_ptr<PData> odata = o.data();
  boost::shared_ptr<PData> outdata = out.data_;
  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset = icount * blocksize_;
    const Complex* cdata = data_->pointer(boffset); 
    const Complex* codata = odata->pointer(boffset); 
    Complex* coutdata = outdata->pointer(boffset);
    zgemm_("N", "N", &l, &n, &m, &one, cdata, &nbasis_, codata, &nbasis_, &zero, coutdata, &nbasis_); 
  }
  out.ndim_ = l;
  out.mdim_ = n; 
  return out;
};


PMatrix1e PMatrix1e::operator%(const PMatrix1e& o) const {
  PMatrix1e out(geom_);
  const int unit = 1;
  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);
  const int l = mdim_;
  const int m = ndim_; 
  assert(ndim_ == o.ndim());
  const int n = o.mdim(); 
  const boost::shared_ptr<PData> odata = o.data();
  boost::shared_ptr<PData> outdata = out.data_;

  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset = icount * blocksize_;
    const Complex* cdata = data_->pointer(boffset);
    const Complex* codata = odata->pointer(boffset); 
    Complex* coutdata = outdata->pointer(boffset);
    zgemm_("C", "N", &l, &n, &m, &one, cdata, &nbasis_, codata, &nbasis_, &zero, coutdata, &nbasis_); 
  }

  out.ndim_ = l;
  out.mdim_ = n;
  return out;
}


void PMatrix1e::diagonalize(double* eig) {
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
    for (int i = 0; i != nbasis_; ++i) {
      for (int j = 0; j != nbasis_; ++j, ++index) {
        cout << setw(14) << setprecision(2) << (*data_)[index] << " ";
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
    for (int i = 0; i != nbasis_; ++i) {
      for (int j = 0; j != nbasis_; ++j, ++index) {
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


