//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_matrix1e_h
#define __src_scf_matrix1e_h

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
    Matrix1e(const std::shared_ptr<Geometry>); 
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
    Matrix1e operator%(const Matrix1e&) const; // caution
    Matrix1e operator^(const Matrix1e&) const; // caution
    Matrix1e operator+(const Matrix1e&) const;
    Matrix1e& operator+=(const Matrix1e&);
    Matrix1e& operator=(const Matrix1e&);
    Matrix1e operator-(const Matrix1e&) const;

    void daxpy(const double, const Matrix1e&);
    void daxpy(const double, const std::shared_ptr<Matrix1e>);
    const double ddot(const Matrix1e&) const;
    const double ddot(const std::shared_ptr<Matrix1e>) const;
    const double rms() const;
    const double trace() const;

    void print(const std::string in = "", const int size = 10) const;
};

#endif

