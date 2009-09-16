//
// Author: Toru Shiozaki
// Date  : April 2009
//
#ifndef __src_rysint_int2d_h
#define __src_rysint_int2d_h 

#include <vector>
#include <map>
#include <src/rysint/scalelist.h>

class Int2D {

  protected:

    /// for recursion
    double C00_[20]; 
    double D00_[20]; 
    double B00_[20]; 
    double B10_[20];
    double B01_[20];

    // two index intermediate
    int rank_;
    double* data_;
    int datasize_;

    ScaleList scale_;

  public:

    Int2D(const double*, const double*, const int, const int, double*, 
          void (*vrrfunc)(double*, const double*, const double*, const double*, const double*, const double*));
    Int2D() {};
    ~Int2D();

    const double* data() const { return data_; };
    const int datasize() const { return datasize_; };

    void scale_data(const double* a, const double c) {
      scale_.scalefunc[rank_](data_, a, c, data_, datasize_);
    };

    void scale_data_t(double* target, const double* a, const double c) {
      scale_.scalefunc[rank_](target, a, c, data_, datasize_);
    };
};


#endif
