//
// Author : std::complex<double>oru Shiozaki
// Date   : July 2009
//

#ifndef __src_util_pmatrix1e_h
#define __src_util_pmatrix1e_h

#define DPI 3.14159265358979323846 

#include <vector>
#include <complex>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <src/pscf/pgeometry.h>
#include <src/pscf/f77.h>
#include <src/util/pdata.h>
#include <boost/shared_ptr.hpp>

class PMatrix1e {
  protected:
    boost::shared_ptr<PGeometry> geom_;
    const int nbasis_;
    // # of row
    int ndim_;
    // # of column
    int mdim_;
    const int blocksize_;
    const int totalsize_;

    virtual void init();
    virtual void computebatch(const std::vector<boost::shared_ptr<Shell> >&,
        const int, const int, const int, const int) { assert(false); };

    boost::shared_ptr<PData> data_;

  public:
    PMatrix1e(const boost::shared_ptr<PGeometry>);
    PMatrix1e(const boost::shared_ptr<PGeometry>, const int ldn, const int ldm);
    // Constructing PMatrix1e while increasing the leading dimension, or ndim_.
    PMatrix1e(const boost::shared_ptr<PMatrix1e> source, const int ldn);
    // Constructing PMatrix1e while reducing the number of columns
    PMatrix1e(const boost::shared_ptr<PMatrix1e> source, const std::pair<int, int> mcut);
    // Constructing Pmatrix merging two matrices
    PMatrix1e(const boost::shared_ptr<PMatrix1e>, const boost::shared_ptr<PMatrix1e>);
    ~PMatrix1e(); 

    PMatrix1e operator*(const PMatrix1e&) const;
    PMatrix1e operator%(const PMatrix1e&) const; // caution
    PMatrix1e operator+(const PMatrix1e&) const;
    PMatrix1e operator-(const PMatrix1e&) const;
    PMatrix1e& operator=(const PMatrix1e&);
    PMatrix1e& operator+=(const PMatrix1e&);

    PMatrix1e ft() const; 
    PMatrix1e bft() const; 

    void set_geom(boost::shared_ptr<PGeometry> a) {geom_ = a; };

    void hermite();
    void real();
    void conj();
    void conj_transpose();
    void scale(const std::complex<double>);
    const std::complex<double>* bp(const int k) const { return data_->pointer((k + geom_->K()) * blocksize_); };
    std::complex<double>* bpw(const int k) { return data_->pointer((k + geom_->K()) * blocksize_); };

    const int nbasis() const { return nbasis_; };
    const int K() const { return geom_->K(); };
    const int L() const { return geom_->L(); };
    const int S() const { return geom_->S(); };
    const double A() const { return geom_->A(); };
    boost::shared_ptr<PData> data() const { return data_; };
    const int ndim() const { return ndim_; };
    const int mdim() const { return mdim_; };
    const int blocksize() const { return blocksize_; };
    const int totalsize() const { return totalsize_; };

    const boost::shared_ptr<PGeometry> geom() const { return geom_; };
    void diagonalize(double*);
    boost::shared_ptr<PMatrix1e> inverse() const;
    void svd(boost::shared_ptr<PMatrix1e> U, boost::shared_ptr<PMatrix1e> V);

    void print() const;
    void rprint(const int precision=6) const;

    void zaxpy(const std::complex<double>, const PMatrix1e&);
    void zaxpy(const std::complex<double>, const boost::shared_ptr<PMatrix1e>);
    const std::complex<double> zdotc(const PMatrix1e&) const;
    const std::complex<double> zdotc(const boost::shared_ptr<PMatrix1e>) const;
   
    const double rms() const;
    const double trace() const;

    std::pair<boost::shared_ptr<PMatrix1e>, boost::shared_ptr<PMatrix1e> >
      split(const int nrow1, const int nrow2) const;
    boost::shared_ptr<PMatrix1e> merge(const boost::shared_ptr<PMatrix1e>) const;

    // AO-to-MO integral transformation...
    boost::shared_ptr<PMatrix1e> mo_transform(boost::shared_ptr<PMatrix1e>,
                             const int istart, const int ifence,
                             const int jstart, const int jfence) const;

}; 

#endif
