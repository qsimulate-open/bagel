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

    double* data_;
    virtual void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int);
    virtual void init();

  public:
    Matrix1e() : nbasis_(0), ndim_(0), mdim_(0) {};
    Matrix1e(const std::shared_ptr<Geometry>); 
    Matrix1e(const std::shared_ptr<Geometry>, const int n, const int m);
    Matrix1e(const Matrix1e&); 
    ~Matrix1e();

    const std::shared_ptr<Geometry> geom() const { return geom_; };

    const int ndim() const { return ndim_; }; 
    const int mdim() const { return mdim_; }; 
    double* data() const { return data_; };
    double& element(int i, int j) { return data_[i+j*ndim_]; };
    double* element_ptr(int i, int j) { return data_+i+j*ndim_; };

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

    void add_diag(const double a, const int i, const int j)
      { for (int ii = i; ii != j; ++ii) data_[ii+ii*nbasis_] += a; };

    void zero() { std::fill(data_, data_+nbasis_*nbasis_, 0.0); };
    void unit() { std::fill(data_, data_+nbasis_*nbasis_, 0.0);
                  for (int i = 0; i != ndim_; ++i) data_[i+i*nbasis_] = 1.0; assert(ndim_ == mdim_);};
    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const Matrix1e& s);

    void print(const std::string in = "", const int size = 10) const;
};

#endif

