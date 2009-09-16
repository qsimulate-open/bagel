//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/matrix1e.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/scf/f77.h>
#include <cassert>
#include <cmath>

using namespace std;

typedef boost::shared_ptr<Geometry> RefGeometry;
typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Shell> RefShell;

Matrix1e::Matrix1e(const RefGeometry geom) : geom_(geom), nbasis_(geom->nbasis()) {
  data_ = new double[nbasis_ * nbasis_]; 
  mdim_ = ndim_ = nbasis_;
  fill(data_, data_ + nbasis_ * nbasis_, 0.0);
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
  const int unit = 1;
  const double one = 1.0; 
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();

  daxpy_(&size, &one, odata, &unit, data_, &unit); 

  return *this; 
}


Matrix1e& Matrix1e::operator=(const Matrix1e& o) {
  const int unit = 1;
  const double one = 1.0; 
  const int size = nbasis_ * nbasis_;
  const double* odata = o.data();

  dcopy_(&size, odata, &unit, data_, &unit); 

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
  const int unit = 1;
  const double one = 1.0;
  const double zero = 0.0;
  const int l = ndim_;
  const int m = mdim_; 
  assert(mdim_ == o.ndim());
  const int n = o.mdim(); 
  const double* odata = o.data();
  double* outdata = out.data_;

  dgemm_("N", "N", &l, &n, &m, &one, data_, &nbasis_, odata, &nbasis_, &zero, outdata, &nbasis_); 

  out.ndim_ = l;
  out.mdim_ = n; 
    
  return out;
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


void Matrix1e::daxpy(const double a, const boost::shared_ptr<Matrix1e> o) {
  const int size = nbasis_ * nbasis_;
  const int unit = 1;
  const double* odata = o->data();
  daxpy_(&size, &a, odata, &unit, data_, &unit); 
}


const double Matrix1e::ddot(const Matrix1e& o) const {
  const int size = nbasis_ * nbasis_;
  const int unit = 1;
  const double* odata = o.data();
  return ddot_(&size, data_, &unit, odata, &unit); 
}


const double Matrix1e::ddot(const boost::shared_ptr<Matrix1e> o) const {
  const int size = nbasis_ * nbasis_;
  const int unit = 1;
  const double* odata = o->data();
  return ddot_(&size, data_, &unit, odata, &unit); 
}


const double Matrix1e::rms() const {
  return ::sqrt(ddot(*this) / (ndim_ * mdim_));
}


const double Matrix1e::trace() const {
  double out = 0.0;
  for (int i = 0; i != ndim_; ++i) out += data_[i * nbasis_ + i]; 
  return out;
}


void Matrix1e::print(const string name, const int size) const { 
 
  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != size; ++i) {
    for (int j = 0; j != size; ++j) {
      cout << fixed << setw(8) << setprecision(5) << data_[j * nbasis_ + i]  << " "; 
    }
    cout << endl;
  }

} 

