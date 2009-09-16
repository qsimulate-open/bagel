//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_util_pdata_h
#define __src_util_pdata_h

#include <complex>
#include <algorithm>

class PData {
  protected:
    std::complex<double>* data_;
    const int length_;

  public:
    PData(int i) : length_(i) { data_ = new std::complex<double>[i];
                                const std::complex<double> zero(0.0, 0.0);
                                std::fill(data_, data_ + i, zero);
                              };
    ~PData() { delete[] data_; };

    std::complex<double>& operator[] (int i) { assert(i < length_ && i >= 0); return data_[i]; };

    std::complex<double>* front() { return data_; };
    std::complex<double>* pointer(size_t j) { return data_ + j; };

};

#endif
