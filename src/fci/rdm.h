//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#ifndef __NEWINT_FCI_RDM_H
#define __NEWINT_FCI_RDM_H

#include <algorithm>
#include <iostream>
#include <iomanip>

class RDM1 {
  protected:
    double* data_;
    const int norb_;
  public:
    RDM1(const int n) : norb_(n) {
      data_ = new double[norb_*norb_];
    };
    RDM1(const RDM1& o) : norb_(o.norb_) {
      data_ = new double[norb_*norb_];
      std::copy(o.data_, o.data_+norb_*norb_, data_);
    }
    ~RDM1() { delete[] data_; };

    double* first() { return data_; };
    double& element(int i, int j) { return *(data_+i+j*norb_); };

    void print() {
      for (int i = 0; i != norb_; ++i) {
        for (int j = 0; j != norb_; ++j) {
          std::cout << std::setw(12) << std::setprecision(7) << element(j,i); 
        }
        std::cout << std::endl;
      }
    };
};

class RDM2 {
  protected:
    double* data_;
    const int norb_;
  public:
    RDM2(const int n) : norb_(n) {
      data_ = new double[norb_*norb_*norb_*norb_];
    };
    RDM2(const RDM2& o) : norb_(o.norb_) {
      data_ = new double[norb_*norb_*norb_*norb_];
      std::copy(o.data_, o.data_+norb_*norb_*norb_*norb_, data_);
    }
    ~RDM2() { delete[] data_; };

    double* first() { return data_; };
    double& element(int i, int j) { return *(data_+i+j*norb_*norb_); };
    double& element(int i, int j, int k, int l) { return *(data_+i+norb_*(j+norb_*(k+norb_*l))); };

    void print() {
      std::cout << "printing 2rdm" << std::endl;
      for (int i = 0; i != norb_; ++i) {
        for (int j = 0; j != norb_; ++j) {
          for (int k = 0; k != norb_; ++k) {
            for (int l = 0; l != norb_; ++l) {
              if (std::abs(element(l,k,j,i)) > 1.0e-2) std::cout << std::setw(3) << l << std::setw(3)
                    << k << std::setw(3) << j << std::setw(3) << i
                    << std::setw(12) << std::setprecision(7) << element(l,k,j,i) << std::endl;
            }
          }
        }
      }
    };
};

#endif
