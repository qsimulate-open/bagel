//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_pscf_diis_h
#define __src_pscf_diis_h

#include <list>
#include <cassert>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <complex>
#include <src/pscf/f77.h>

// boost::shared_ptr<T> is assumed to be a shared_pointer of some class
// which have zaxpy and zdotc functions.
// T must have geom() function that returns Geometry. 

// Complex version

template <class T>
class PDIIS {
  typedef boost::shared_ptr<T> RefT;
  typedef std::list<std::pair<RefT, RefT> > Container_type_;
  typedef typename Container_type_::iterator iterator;
  typedef std::complex<double> Complex;

  protected:
    const int ndiis_;
    const int nld_;
    const double dumping_;

    Container_type_ data_;

    Complex* matrix_;
    Complex* matrix_save_;
    Complex* coeff_;
    Complex* work_;
    int* ipiv_;
    const int lwork_;
  
  public:
    PDIIS(const int ndiis, const double dump) : ndiis_(ndiis), nld_(ndiis + 1), lwork_(nld_ * nld_), dumping_(dump) {
      matrix_ = new Complex[nld_ * nld_]; 
      matrix_save_ = new Complex[nld_ * nld_]; 
      coeff_ = new Complex[nld_]; 
      work_ = new Complex[lwork_];
      ipiv_ = new int[nld_];
    };  

    ~PDIIS() {
      delete[] matrix_;
      delete[] matrix_save_;
      delete[] coeff_;
      delete[] work_;
      delete[] ipiv_;
    };

    RefT extrapolate(const std::pair<RefT, RefT> input) {
      RefT v = input.first;
      RefT e = input.second;
      data_.push_back(input);
      const int unit = 1;
      const int size = nld_ * nld_;

      if (data_.size() > ndiis_) {
        data_.pop_front();
        for (int i = 1; i != ndiis_; ++i) { 
          for (int j = 1; j != ndiis_; ++j)
            matrix_[(j - 1) + (i - 1) * nld_] = matrix_save_[j + i * nld_]; 
        }
      } else if (data_.size() != 1) {
        zcopy_(&size, matrix_save_, &unit, matrix_, &unit);
      } 
      const int cnum = data_.size();
      iterator data_iter = data_.begin();

      for (int i = 0; i != cnum - 1; ++i, ++data_iter) 
        matrix_[(cnum - 1) + i * nld_] = matrix_[i + (cnum - 1) * nld_] = (e->zdotc(*(data_iter->second))).real(); 
      matrix_[(cnum - 1) + (cnum - 1) * nld_] = (e->zdotc(e)).real();
      for (int i = 0; i != cnum; ++i)
        matrix_[cnum + i * nld_] = matrix_[i + cnum * nld_] = -1.0;
      matrix_[cnum + cnum * nld_] = 0.0;
      for (int i = 0; i != cnum; ++i) coeff_[i] = 0.0;
      coeff_[cnum] = -1.0;

      zcopy_(&size, matrix_, &unit, matrix_save_, &unit); 

      const int cdim = cnum + 1;
      int info_in_diis;
      zhesv_("U", &cdim, &unit, matrix_, &nld_, ipiv_, coeff_, &nld_, work_, &lwork_, &info_in_diis);
      if(info_in_diis != 0) throw std::runtime_error("DHESV in DIIS failed");

      RefT out(new T(e->geom())); 

      data_iter = data_.begin();
      for (int i = 0; i != cnum; ++i, ++data_iter)
        out->zaxpy((coeff_[i]).real() * (1.0 - dumping_), *(data_iter->first));
      out->zaxpy(dumping_, *(input.first));
      out->real();

      return out;
    };

};

#endif

