//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_matrix1e_h
#define __src_scf_matrix1e_h

#include <cassert>
#include <src/scf/shell.h>
#include <src/scf/geometry.h>
#include <string>
#include <algorithm>
#include <memory>

class Matrix1e {
  protected:
    std::shared_ptr<Geometry> geom_;
    int nbasis_;
    int ndim_;
    int mdim_;

    std::unique_ptr<double[]> data_;
    virtual void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int);
    virtual void init();

  public:
    Matrix1e() : nbasis_(0), ndim_(0), mdim_(0) {};
    Matrix1e(const std::shared_ptr<Geometry>); 
    Matrix1e(const std::shared_ptr<Geometry>, const int n, const int m);
    Matrix1e(const Matrix1e&); 
    ~Matrix1e();

    const std::shared_ptr<Geometry> geom() const { return geom_; };

    const int nbasis() const { return nbasis_; };
    const int size() const { return nbasis_*nbasis_; };
    const int ndim() const { return ndim_; }; 
    const int mdim() const { return mdim_; }; 
    double* data() const { return data_.get(); };
    double& data(const size_t i) { return *(data_.get()+i); };
    double& element(int i, int j) { return *element_ptr(i, j); };
    double* element_ptr(int i, int j) { return data()+i+j*ndim_; };

    void symmetrize();
    void diagonalize(double*);
    void inverse();

    Matrix1e operator*(const Matrix1e&) const;
    Matrix1e& operator*=(const Matrix1e&);
    Matrix1e operator*(const double& a) const;
    Matrix1e& operator*=(const double& a);
    Matrix1e operator%(const Matrix1e&) const; // caution
    Matrix1e operator^(const Matrix1e&) const; // caution
    Matrix1e operator+(const Matrix1e&) const;
    Matrix1e& operator+=(const Matrix1e&);
    Matrix1e& operator-=(const Matrix1e&);
    Matrix1e& operator=(const Matrix1e&);
    Matrix1e operator-(const Matrix1e&) const;

    std::shared_ptr<Matrix1e> clone() const {
      std::shared_ptr<Matrix1e> out(new Matrix1e(geom_, ndim_, mdim_));
      return out;
    };
    // returns exp(*this)
    std::shared_ptr<Matrix1e> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix1e> log(const int deg = 6) const;

    void daxpy(const double, const Matrix1e&);
    void daxpy(const double, const std::shared_ptr<Matrix1e>);
    const double ddot(const Matrix1e&) const;
    const double norm() const { return std::sqrt(ddot(*this)); };
    const double ddot(const std::shared_ptr<Matrix1e>) const;
    const double rms() const;
    const double trace() const;
    void dscal(const double a) { dscal_(size(), a, data(), 1); };

    void add_diag(const double a, const int i, const int j)
      { for (int ii = i; ii != j; ++ii) element(ii,ii) += a; };

    void zero() { std::fill(data(), data()+nbasis_*nbasis_, 0.0); };
    void unit() { std::fill(data(), data()+nbasis_*nbasis_, 0.0);
                  for (int i = 0; i != ndim_; ++i) element(i,i) = 1.0; assert(ndim_ == mdim_);};
    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const Matrix1e& s);

    void print(const std::string in = "", const int size = 10) const;

    double orthog(const std::list<std::shared_ptr<Matrix1e> > o) {
      for (auto iter = o.begin(); iter != o.end(); ++iter) daxpy(-ddot(*iter), *iter);
      const double n = norm(); *this *= 1.0/n;
      return n;
    };
};

#endif

