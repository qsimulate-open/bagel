//
// BAGEL - Parallel electron correlation program.
// Filename: pmatrix1e.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <complex>
#include <src/util/constants.h>
#include <src/pscf/pmatrix1e.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;

typedef std::complex<double> Complex;


// I have supposed that ndim_ is the leading dimension in the data.

PMatrix1e::PMatrix1e(const shared_ptr<PGeometry> pg)
: geom_(pg), nbasis_(pg->nbasis()),
  blocksize_(pg->nbasis() * pg->nbasis()), totalsize_((2 * pg->K() + 1) * pg->nbasis() * pg->nbasis()) {

  mdim_ = ndim_ = nbasis_;
  shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
}


PMatrix1e::PMatrix1e(const shared_ptr<PGeometry> pg, const int ldn, const int ldm)
: geom_(pg), nbasis_(pg->nbasis()),
  blocksize_(ldn * ldm), totalsize_((2 * pg->K() + 1) * ldn * ldm) {

  mdim_ = ldm;
  ndim_ = ldn;
  shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
}


// Changes the leading dimension, filling zero if expanding.
PMatrix1e::PMatrix1e(const shared_ptr<PMatrix1e> source, const int ldn)
: geom_(source->geom()), nbasis_(source->nbasis()),
  blocksize_(ldn * source->mdim()), totalsize_((2 * source->K() + 1) * ldn * source->mdim()) {

  assert(ldn >= nbasis_);

  ndim_ = ldn;
  mdim_ = source->mdim();
  shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
  const shared_ptr<PData> sdata = source->data();
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


// Changes number of columns; keeps only [mcut.first mcut.second) columns.
PMatrix1e::PMatrix1e(const shared_ptr<PMatrix1e> source, pair<int, int> mcut)
: geom_(source->geom()), nbasis_(source->nbasis()),
  blocksize_(source->ndim() * (mcut.second-mcut.first)),
  totalsize_((2 * source->K() + 1) * source->ndim() * (mcut.second-mcut.first)) {

  assert(mcut.second <= source->mdim());
  assert(mcut.first >= 0 && mcut.first < mcut.second);

  ndim_ = source->ndim();
  mdim_ = mcut.second-mcut.first;
  shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
  const shared_ptr<PData> sdata = source->data();
  const int sblocksize = source->blocksize();

  for (int i = -K(); i != max(K(), 1); ++i) {
    complex<double>* current = data_->front() + (i + K()) * blocksize_;
    complex<double>* csource = sdata->front() + (i + K()) * sblocksize;

    copy(csource+ndim_*mcut.first, csource+ndim_*mcut.second, current);
  }
}


PMatrix1e::PMatrix1e(const shared_ptr<PMatrix1e> source1, const shared_ptr<PMatrix1e> source2)
: geom_(source1->geom()), nbasis_(source1->nbasis()), ndim_(source1->ndim()), mdim_(source1->mdim()+source2->mdim()),
  blocksize_(source1->ndim() * (source1->mdim()+source2->mdim())),
  totalsize_((2*source1->K()+1)*source1->ndim()*(source1->mdim()+source2->mdim())) {

  if (source1->ndim() != source2->ndim()) {
    cout << " source 1 dimension: " << source1->ndim() << " " << source1->mdim() << " " << source1->nbasis() << " " << source1->blocksize() << endl;
    cout << " source 2 dimension: " << source2->ndim() << " " << source2->mdim() << " " << source2->nbasis() << " " << source2->blocksize() << endl;
    throw runtime_error("source 1 and source 2 should have the same ndim");
  }

  shared_ptr<PData> tmp(new PData(totalsize_));
  data_ = tmp;
  const shared_ptr<PData> sdata1 = source1->data();
  const shared_ptr<PData> sdata2 = source2->data();
  const size_t sblocksize1 = source1->blocksize();
  const size_t sblocksize2 = source2->blocksize();

  for (int i = -K(); i != max(K(), 1); ++i) {
    complex<double>* current = data_->front() + (i + K()) * blocksize_;
    complex<double>* csource1 = sdata1->front() + (i + K()) * sblocksize1;
    complex<double>* csource2 = sdata2->front() + (i + K()) * sblocksize2;

    copy(csource1, csource1+sblocksize1, current);
    copy(csource2, csource2+sblocksize2, current+sblocksize1);
  }
}


PMatrix1e::~PMatrix1e() {
}


PMatrix1e PMatrix1e::operator+(const PMatrix1e& add) const {
  assert(add.totalsize() == totalsize_);
  PMatrix1e out(geom_, ndim_, mdim_);
  for (int j = 0; j != totalsize_; ++j) (*out.data_)[j] = (*add.data_)[j] + (*data_)[j];
  return out;
}


PMatrix1e PMatrix1e::operator-(const PMatrix1e& sub) const {
  assert(sub.totalsize() == totalsize_);
  PMatrix1e out(geom_, ndim_, mdim_);
//#pragma omp parallel for
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
  zaxpy_(&totalsize_, &one, source.data_->cfront(), &unit, data_->front(), &unit);
  return *this;
}


PMatrix1e PMatrix1e::ft() const {
  PMatrix1e out(geom_);

  const int unit = 1;
  for (int k = -K(); k < max(K(), 1); ++k) {
    const Complex kq(0.0, K() != 0 ? (pi__ * k) / K() : 0.0);
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
        const Complex kq(0.0, (pi__ * k) / K());
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
  typedef shared_ptr<const Atom> RefAtom;
  typedef shared_ptr<const Shell> RefShell;

  const std::vector<RefAtom> atoms = geom_->atoms();

  const std::vector<std::vector<int>> offsets = geom_->offsets();
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
        array<double,3> disp;
        disp[0] = disp[1] = 0.0;
        disp[2] = geom_->A() * m;
        const RefAtom catom1(new Atom(*(atoms[iatom1]), disp));
        const int numshell1 = catom1->shells().size();
        const std::vector<int> coffset1 = offsets[iatom1];
        const std::vector<RefShell> shell1 = catom1->shells();

        for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
          const int offset0 = coffset0[ibatch0];
          const shared_ptr<const Shell> b0 = shell0[ibatch0];
          for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
            const int offset1 = coffset1[ibatch1];
            const shared_ptr<const Shell> b1 = shell1[ibatch1];
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
  const shared_ptr<PData> odata = o.data();
  shared_ptr<PData> outdata = out.data();

  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset1 = icount * blocksize_;
    const int boffset2 = icount * o.blocksize();
    const int boffset3 = icount * out.blocksize();
    const Complex* cdata = data_->pointer(boffset1);
    const Complex* codata = odata->pointer(boffset2);
    Complex* coutdata = outdata->pointer(boffset3);
    zgemm3m_("N", "N", &l, &n, &m, &one, cdata, &ndim_, codata, &mdim_, &zero, coutdata, &ndim_);
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
  const shared_ptr<PData> odata = o.data();

  PMatrix1e out(geom_, l, n);
  shared_ptr<PData> outdata = out.data();

  for (int i = -K(); i <= K(); ++i) {
    const int icount = i + K();
    const int boffset1 = icount * blocksize_;
    const int boffset2 = icount * o.blocksize();
    const int boffset3 = icount * out.blocksize();
    const Complex* cdata = data_->pointer(boffset1);
    const Complex* codata = odata->pointer(boffset2);
    Complex* coutdata = outdata->pointer(boffset3);
    zgemm3m_("C", "N", &l, &n, &m, &one, cdata, &ndim_, codata, &ndim_, &zero, coutdata, &mdim_);
  }

  return out;
}


tuple<shared_ptr<PMatrix1e>, shared_ptr<PMatrix1e>> PMatrix1e::svd() {
  auto U = make_shared<PMatrix1e>(geom_, ndim_, ndim_);
  auto V = make_shared<PMatrix1e>(geom_, mdim_, mdim_);

  const int lwork = 5 * min(ndim_, mdim_)*min(ndim_, mdim_) + 5*max(ndim_, mdim_);
  complex<double>* work = new complex<double>[lwork];
  double* rwork = new double[lwork];
  double* S = new double[min(ndim_, mdim_)];
  int* iwork = new int[min(ndim_, mdim_)*8];
/*
  SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, RWORK, INFO )
 */
  for (int i = -K(); i < max(K(), 1); ++i) {
    complex<double>* cblock = data_->pointer(blocksize_ * (i + K()));
    complex<double>* ublock = U->data()->pointer(U->blocksize() * (i + K()));
    complex<double>* vblock = V->data()->pointer(V->blocksize() * (i + K()));
    int info = 0;
#define USE_DC
#ifdef USE_DC
    zgesdd_("A", &ndim_, &mdim_, cblock, &ndim_, S, ublock, &ndim_, vblock, &mdim_, work, &lwork, rwork, iwork, &info);
#else
    zgesvd_("A", "A", &ndim_, &mdim_, cblock, &ndim_, S, ublock, &ndim_, vblock, &mdim_, work, &lwork, rwork, &info);
#endif
    if (info != 0) throw runtime_error("zgesvd failed in PMatrix1e::svd");
#if 0
    for (int j = 0; j != min(ndim_, mdim_); ++j) cout << setprecision(10) << fixed << S[j] << " ";
    cout << endl;
#endif
  }
  delete[] S;
  delete[] work;
  delete[] rwork;
  delete[] iwork;

  return make_tuple(U,V);
}


void PMatrix1e::diagonalize(double* eig) {
  assert(ndim_ == mdim_);
  const int lwork = ndim_ * 6;

//#pragma omp parallel for
  for (int i = -K(); i <= K(); ++i) {
    Complex* work = new Complex[lwork];
    double* rwork = new double[ndim_ * 3 - 2];
    const int icount = i + K();
    const int boffset = icount * blocksize_;
    const int eoffset = icount * ndim_;
    Complex* cdata = data_->pointer(boffset);
    double* ceig = &eig[eoffset];

    int info_diagonalize = 0;
    zheev_("V", "L", &ndim_, cdata, &ndim_, ceig, work, &lwork, rwork, &info_diagonalize);
    assert(info_diagonalize == 0);
    delete[] work;
    delete[] rwork;
  }

}


std::shared_ptr<PMatrix1e> PMatrix1e::inverse() const {
//        SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  assert(ndim_ == mdim_);
  PMatrix1e cpy(geom_, ndim_, mdim_);
  cpy = *this;
  shared_ptr<PMatrix1e> out(new PMatrix1e(geom_, ndim_, mdim_));

//#pragma omp parallel for
  for (int i = -K(); i < std::max(K(), 1); ++i) {
    int* ipiv = new int[ndim_];

    const int icount = i + K();
    const int boffset = icount * blocksize_;
    Complex* cdata = cpy.data()->pointer(boffset);

    Complex* odata = out->data()->pointer(boffset);
    // setting unit matrix.
    for (int i = 0; i != ndim_; ++i) {
      odata[i + i * ndim_] = 1.0;
    }
    int info;
    zgesv_(&ndim_, &ndim_, cdata, &ndim_, ipiv, odata, &ndim_, &info);

    if (info != 0) throw runtime_error("PMatrix1e::inverse() failed.");

    delete[] ipiv;
  }
  return out;
}


void PMatrix1e::print() const {
  int index = 0;
  for (int k = -K(); k <= K(); ++k) {
    cout << "K = " << setw(2) << k << ":" << endl;
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j, ++index) {
        cout << setw(14) << setprecision(2) << (*data_)[index] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}


void PMatrix1e::rprint(const int precision) const {
  int index = 0;
  for (int k = -K(); k <= K(); ++k) {
    cout << "K = " << setw(2) << k << ":" << endl;
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j, ++index) {
        cout << setw(precision+3) << setprecision(precision) << (*data_)[index].real() << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}


void PMatrix1e::hermite() {
  assert(ndim_ == mdim_);
//#pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    const int kcount = k + K();
    const int koffset = kcount * blocksize_;
    Complex* dat = data_->pointer(koffset);
    for (int i = 0; i != ndim_; ++i) {
      for (int j = 0; j != ndim_; ++j) {
        const Complex a = dat[j + i * ndim_];
        const Complex b = dat[i + j * ndim_];
        dat[j + i * ndim_] = (a + std::conj(b)) * 0.5;
        dat[i + j * ndim_] = (b + std::conj(a)) * 0.5;
      }
    }
  }
}


void PMatrix1e::real() {
//#pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    const int kcount = k + K();
    const int koffset = kcount * blocksize_;
    Complex* dat = data_->pointer(koffset);
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j) {
        dat[j + i * ndim_] = (dat[j + i * ndim_]).real();
      }
    }
  }
}


void PMatrix1e::conj() {
//#pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    const int kcount = k + K();
    const int koffset = kcount * blocksize_;
    Complex* dat = data_->pointer(koffset);
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j) {
        dat[j + i * ndim_] = std::conj(dat[j + i * ndim_]);
      }
    }
  }
}


void PMatrix1e::conj_transpose() {
  assert(ndim_ == mdim_);
//#pragma omp parallel for
  for (int k = -K(); k <= K(); ++k) {
    Complex* work = new Complex[blocksize_];
    const int kcount = k + K();
    const int koffset = kcount * blocksize_;
    Complex* dat = data_->pointer(koffset);
    copy(dat, dat+blocksize_, work);
    for (int i = 0; i != mdim_; ++i) {
      for (int j = 0; j != ndim_; ++j) {
        dat[j + i * ndim_] = std::conj(work[i + j * ndim_]);
      }
    }
    delete[] work;
  }
}


void PMatrix1e::scale(const Complex a) {
  const int unit = 1;
  zscal_(&totalsize_, &a, data_->front(), &unit);
}


void PMatrix1e::ax_plus_y(const Complex a, const PMatrix1e& o) {
  const int unit = 1;
  const Complex* odata = o.data()->front();
  zaxpy_(&totalsize_, &a, odata, &unit, data_->front(), &unit);
}


void PMatrix1e::ax_plus_y(const Complex a, const shared_ptr<PMatrix1e> o) {
  const int unit = 1;
  const Complex* odata = o->data()->front();
  zaxpy_(&totalsize_, &a, odata, &unit, data_->front(), &unit);
}


const Complex PMatrix1e::dot_product(const PMatrix1e& o) const {
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


const Complex PMatrix1e::dot_product(const shared_ptr<PMatrix1e> o) const {
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


double PMatrix1e::rms() const {
  const int unit = 1;
  Complex zdot;
#ifndef ZDOT_RETURN
  zdotc_(&zdot, &totalsize_, data_->front(), &unit, data_->front(), &unit);
#else
  zdot = zdotc_(&totalsize_, data_->front(), &unit, data_->front(), &unit);
#endif
  return ::sqrt(abs(zdot) / totalsize_);
}


double PMatrix1e::trace() const {
  assert(ndim_ == mdim_);
  double out = 0.0;
  int scount = 0;
  for (int s = -L(); s <= L(); ++s, ++scount) {
    const int boffset =  scount * blocksize_;
    for (int i = 0; i != ndim_; ++i)
      out += ((*data_)[boffset + i * ndim_ + i]).real();
  }
  return out;
}


pair<shared_ptr<PMatrix1e>, shared_ptr<PMatrix1e>> PMatrix1e::split(const int nrow1, const int nrow2) const {
  shared_ptr<PMatrix1e> out1(new PMatrix1e(geom_, nrow1, mdim_));
  shared_ptr<PMatrix1e> out2(new PMatrix1e(geom_, nrow2, mdim_));

  assert(nrow1+nrow2 == ndim_);
  assert(blocksize_ == out1->blocksize() + out2->blocksize());
  assert(blocksize_ == ndim_ * mdim_);

  const Complex* source = data_->cfront();
  Complex* data1 = out1->data()->front();
  Complex* data2 = out2->data()->front();

  for (int i = -K(); i <= K(); ++i) {
    for (int m = 0; m != mdim_; ++m, data1+=out1->ndim(), data2+=out2->ndim(), source+=ndim_) {
      copy(source,       source+nrow1,       data1);
      copy(source+nrow1, source+nrow1+nrow2, data2);
    }
  }

  return make_pair(out1, out2);
}


shared_ptr<PMatrix1e> PMatrix1e::merge(const shared_ptr<PMatrix1e> in2) const {
  assert(mdim_ == in2->mdim_);
  shared_ptr<PMatrix1e> out(new PMatrix1e(geom_, ndim_+in2->ndim_, mdim_));

  Complex* outdata = out->data_->front();
  const Complex* data1 = data_->cfront();
  const Complex* data2 = in2->data_->cfront();

  for (int i = -K(); i <= K(); ++i) {
    for (int m = 0; m != mdim_; ++m, data1+=ndim_, data2+=in2->ndim_, outdata+=out->ndim_) {
      copy(data1, data1+ndim_,      outdata);
      copy(data2, data2+in2->ndim_, outdata+ndim_);
    }
  }

  return out;
}


shared_ptr<PMatrix1e> PMatrix1e::mo_transform(shared_ptr<PMatrix1e> coeff,
                         const int istart, const int ifence,
                         const int jstart, const int jfence) const {
  assert(coeff->ndim() == ndim_);
  assert(coeff->ndim() == mdim_);

  const int isize = ifence-istart;
  const int jsize = jfence-jstart;

  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);

  shared_ptr<PMatrix1e> out(new PMatrix1e(geom_, jfence-jstart, ifence-istart));

  Complex* tmp = new Complex[ndim_ * isize];
  Complex* tmp2 = new Complex[ndim_ * isize];
  const int kk = std::max(geom_->K(), 1);

  for (int k = -K(); k != kk; ++k) {
    const Complex* data = this->bp(k);
    const Complex* co = coeff->bp(k);
    Complex* outdata = out->bpw(k);

    int iall = 0;
    for (int i = istart*ndim_; i != ifence*ndim_; ++i, ++iall) tmp[iall] = std::conj(co[i]);

    zgemm3m_("N", "N", &ndim_, &isize, &ndim_, &one, data, &ndim_, tmp, &ndim_, &zero, tmp2, &ndim_);

    const int unit = 1;
    for (int i = 0; i != isize; ++i) {
      zgemm3m_("C", "N", &unit, &jsize, &ndim_, &one, tmp2+i*ndim_, &unit, co+jstart*ndim_, &ndim_, &zero, outdata+i*jsize, &unit);
    }
  }
  delete[] tmp;
  delete[] tmp2;
  return out;

}


