//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/matrix1e.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <cassert>
#include <cmath>

using namespace std;

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

Matrix1e::Matrix1e(const RefGeometry geom) : geom_(geom), nbasis_(geom->nbasis()) {
  data_ = new double[nbasis_ * nbasis_]; 
  mdim_ = ndim_ = nbasis_;
  fill(data_, data_ + nbasis_ * nbasis_, 0.0);
}


Matrix1e::Matrix1e(const RefGeometry geom, const int n, const int m) : geom_(geom), nbasis_(geom->nbasis()) {
  ndim_ = n;
  mdim_ = m;
  data_ = new double[nbasis_*nbasis_]; 
  fill(data_, data_ + nbasis_*nbasis_, 0.0);
}


Matrix1e::Matrix1e(const Matrix1e& o) : geom_(o.geom_), nbasis_(o.nbasis_), ndim_(o.ndim_), mdim_(o.mdim_) {
  data_ = new double[nbasis_ * nbasis_]; 
  copy(o.data_, o.data_ + nbasis_*nbasis_, data_);
}


void Matrix1e::init() {

  const vector<RefAtom> atoms = geom_->atoms(); 
  vector<Atom>::const_iterator aiter0, aiter1;
  vector<RefShell>::const_iterator biter0, biter1;

  const vector<vector<int> > offsets = geom_->offsets();

  const int nbasis = geom_->nbasis();

  // only lower half will be stored
  for (int iatom0 = 0; iatom0 != geom_->natom(); ++iatom0) {
    // iatom1 = iatom1;
    const RefAtom catom0 = atoms[iatom0];
    const int numshell0 = catom0->shells().size();
    const vector<int> coffset0 = offsets[iatom0];
    const vector<RefShell> shell0 = catom0->shells();

    for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
      const int offset0 = coffset0[ibatch0]; 
      RefShell b0 = shell0[ibatch0];
      for (int ibatch1 = ibatch0; ibatch1 != numshell0; ++ibatch1) {
        const int offset1 = coffset0[ibatch1]; 
        RefShell b1 = shell0[ibatch1];
        vector<RefShell> input;
        input.push_back(b1);
        input.push_back(b0);

        computebatch(input, offset0, offset1, nbasis);
      }
    } 

    for (int iatom1 = iatom0 + 1; iatom1 != geom_->natom(); ++iatom1) {
      const RefAtom catom1 = atoms[iatom1];
      const int numshell1 = catom1->shells().size();
      const vector<int> coffset1 = offsets[iatom1];
      const vector<RefShell> shell1 = catom1->shells();

      for (int ibatch0 = 0; ibatch0 != numshell0; ++ibatch0) {
        const int offset0 = coffset0[ibatch0]; 
        RefShell b0 = shell0[ibatch0];
        for (int ibatch1 = 0; ibatch1 != numshell1; ++ibatch1) {
          const int offset1 = coffset1[ibatch1]; 
          RefShell b1 = shell1[ibatch1];
          vector<RefShell> input;
          input.push_back(b1);
          input.push_back(b0);

          computebatch(input, offset0, offset1, nbasis);

        }
      } 
    }
  } 

}


Matrix1e::~Matrix1e() {
  delete[] data_;

}


void Matrix1e::computebatch(const std::vector<RefShell>&, const int, const int, const int) {
  assert(false);
}


void Matrix1e::symmetrize() {
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = i + 1; j != nbasis_; ++j) {
      data_[i + j * nbasis_] = data_[j + i * nbasis_]; 
    }
  } 
}


Matrix1e Matrix1e::operator+(const Matrix1e& o) const {
  Matrix1e out(geom_);
  out.ndim_ = ndim_;
  out.mdim_ = mdim_;
  const int unit = 1;
  const double one = 1.0; 
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();
  double* outdata = out.data_;

  daxpy_(&size, &one, data_, &unit, outdata, &unit); 
  daxpy_(&size, &one, odata, &unit, outdata, &unit); 

  return out; 
}


Matrix1e& Matrix1e::operator+=(const Matrix1e& o) {
  daxpy_(nbasis_*nbasis_, 1.0, o.data(), 1, data_, 1); 
  return *this; 
}


Matrix1e& Matrix1e::operator=(const Matrix1e& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  dcopy_(nbasis_*nbasis_, o.data(), 1, data_, 1); 
  return *this; 
}


Matrix1e Matrix1e::operator-(const Matrix1e& o) const {
  Matrix1e out(geom_);
  out.ndim_ = ndim_;
  out.mdim_ = mdim_;
  const int unit = 1;
  const double one = 1.0; 
  const double mone = -1.0; 
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();
  double* outdata = out.data_;

  daxpy_(&size, &one, data_, &unit, outdata, &unit); 
  daxpy_(&size, &mone, odata, &unit, outdata, &unit); 

  return out; 
}


Matrix1e Matrix1e::operator*(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  double* outdata = out.data_;

  dgemm_("N", "N", l, n, m, 1.0, data_, nbasis_, odata, nbasis_, 0.0, outdata, nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n; 
  return out;
}


Matrix1e& Matrix1e::operator*=(const Matrix1e& o) {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  dgemm_("N", "N", l, n, m, 1.0, data_, nbasis_, odata, nbasis_, 0.0, out.data(), nbasis_); 
  *this = out;
}


Matrix1e Matrix1e::operator*(const double& a) const {
  Matrix1e out(*this);
  dscal_(nbasis_*nbasis_, a, out.data(), 1);
  return out;
}

Matrix1e& Matrix1e::operator*=(const double& a) {
  dscal_(nbasis_*nbasis_, a, data_, 1);
  return *this;
}


Matrix1e Matrix1e::operator%(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int unit = 1;
  const double one = 1.0;
  const double zero = 0.0;
  const int l = mdim_;
  const int m = ndim_; 
  assert(ndim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  double* outdata = out.data_;

  dgemm_("T", "N", &l, &n, &m, &one, data_, &nbasis_, odata, &nbasis_, &zero, outdata, &nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n;
    
  return out;
}


Matrix1e Matrix1e::operator^(const Matrix1e& o) const {
  Matrix1e out(geom_);
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.mdim());
  const int n = o.ndim(); 
  const double* odata = o.data();
  double* outdata = out.data_;

  dgemm_("N", "T", l, n, m, 1.0, data_, nbasis_, odata, nbasis_, 0.0, outdata, nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n;
    
  return out;
}


void Matrix1e::diagonalize(double* eig) {
  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
  
  int info_diagonalize = 0;
  const int lwork = nbasis_ * 6;
  double* work = new double[lwork];
  dsyev_("V", "L", &ndim_, data_, &nbasis_, eig, work, &lwork, &info_diagonalize); 

  assert(info_diagonalize == 0);

  delete[] work;
}


void Matrix1e::daxpy(const double a, const Matrix1e& o) {
  const int size = nbasis_ * nbasis_;
  const int unit = 1;
  const double* odata = o.data();
  daxpy_(&size, &a, odata, &unit, data_, &unit); 
}


void Matrix1e::daxpy(const double a, const std::shared_ptr<Matrix1e> o) {
  const int size = nbasis_ * nbasis_;
  const int unit = 1;
  const double* odata = o->data();
  daxpy_(&size, &a, odata, &unit, data_, &unit); 
}


const double Matrix1e::ddot(const Matrix1e& o) const {
  return ddot_(nbasis_*nbasis_, data_, 1, o.data(), 1); 
}


const double Matrix1e::ddot(const std::shared_ptr<Matrix1e> o) const {
  return ddot_(nbasis_*nbasis_, data_, 1, o->data(), 1); 
}


const double Matrix1e::rms() const {
  return ::sqrt(ddot(*this) / (ndim_ * mdim_));
}


const double Matrix1e::trace() const {
  double out = 0.0;
  for (int i = 0; i != ndim_; ++i) out += data_[i * nbasis_ + i]; 
  return out;
}


shared_ptr<Matrix1e> Matrix1e::exp(const int deg) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_, ndim_, mdim_));
  Matrix1e buf(*this);
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = 1.0/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0; 
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim_; ++j) out->element(j,j) += 1.0; 
  return out;
}


shared_ptr<Matrix1e> Matrix1e::log(const int deg) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_, ndim_, mdim_));
  Matrix1e buf(*this);
  for (int j = 0; j != ndim_; ++j) buf.element(j,j) -= 1.0;
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = -static_cast<double>(i-1)/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0; 
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


void Matrix1e::purify_unitary() {
  assert(ndim_ == mdim_);
  Matrix1e buf(*this ^ *this);
  const int lwork = 5*ndim_;
  double* work = new double[lwork];
  double* vec = new double[ndim_];
  int info;
  dsyev_("V", "L", &ndim_, buf.data(), &ndim_, vec, work, &lwork, &info); 
  if (info) throw runtime_error("dsyev failed in Matrix1e::purify_unitary");
  if (vec[0] < 0.95)        cout << "   --- smallest eigenvalue in purify_unitary() " << vec[0] << endl;
  if (vec[ndim_-1] > 1.05)  cout << "   --- largest eigenvalue in purify_unitary() " << vec[ndim_-1] << endl;
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      buf.element(j,i) /= sqrt(sqrt(vec[i]));
    }
  }
  *this *= buf ^ buf;

  // just checking...
  assert(std::abs((*this^*this).norm()-sqrt(static_cast<double>(ndim_))) < 1.0e-10);

  delete[] work;
  delete[] vec;
}


// in-place matrix inverse (practically we use buffer area)
void Matrix1e::inverse() {
  assert(ndim_ == mdim_);
  shared_ptr<Matrix1e> buf = this->clone();
  buf->unit();

  const int lwork = 3*ndim_;
  int info;
  int* ipiv = new int[ndim_];
#if 0
  double* work = new double[lwork];
  dsysv_("U", &ndim_, &ndim_, data_, &ndim_, ipiv, buf->data(), &ndim_, work, &lwork, &info); 
  delete[] work;
#else
  dgesv_(&ndim_, &ndim_, data_, &ndim_, ipiv, buf->data(), &ndim_, &info);
#endif
  if (info) throw runtime_error("dsysv failed in Matrix1e::inverse()");

  copy(buf->data_, buf->data_+nbasis_*nbasis_, data_);

  delete[] ipiv;
}


void Matrix1e::print(const string name, const int size) const { 
 
  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != size; ++i) {
    for (int j = 0; j != size; ++j) {
      cout << fixed << setw(12) << setprecision(8) << data_[j * nbasis_ + i]  << " "; 
    }
    cout << endl;
  }

} 

