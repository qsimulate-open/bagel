//
// Author : Toru Shiozaki
// Date   : September 2009
//

#pragma once
#ifndef __src_util_pmodiagfile_h
#define __src_util_pmodiagfile_h

#define GEMINAL_EXP 1.5

#include <src/util/pfile.h>
#include <src/pscf/f77.h>

template<class T>
class PMODiagFile : public PFile<T> {
  protected:
    const int istart_;
    const int jstart_;
    const int ifence_;
    const int jfence_;

  public:
    PMODiagFile(const long, const int, const int, const int, const int, const int, const bool late_init = false); 
    ~PMODiagFile(); 

    PMODiagFile<T> operator+(const PMODiagFile<T>&) const;
    PMODiagFile<T> operator-(const PMODiagFile<T>&) const;

    void print() const;
    void rprint() const;

    T get_energy_one_amp() const;

};


template<class T>
PMODiagFile<T>::PMODiagFile(const long fsize, const int k,
                            const int istrt, const int ifen,
                            const int jstrt, const int jfen,
                            const bool late_init)
 : PFile<T>(fsize, k, late_init), istart_(istrt), ifence_(ifen), jstart_(jstrt), jfence_(jfen) {


};


template<class T>
PMODiagFile<T>::~PMODiagFile() {
};


template<class T>
PMODiagFile<T> PMODiagFile<T>::operator+(const PMODiagFile<T>& other) const {
  const int k = this->K_;
  const size_t fsize = this->filesize_;
  const int kk4 = std::max(k * k * 4, 1);
  const size_t blocksize = fsize / kk4;
  T* buffer1 = new T[blocksize];
  T* buffer2 = new T[blocksize];
  assert(fsize % kk4 == 0);

  PMODiagFile out(fsize, k, istart_, ifence_, jstart_, jfence_, false);

  size_t current = 0lu;
  for (int iblock = 0; iblock != kk4; ++iblock, current += blocksize) {
    this->get_block(current, blocksize, buffer1); 
    other.get_block(current, blocksize, buffer2);

    #pragma omp parallel for schedule(dynamic, 128)
    for (int i = 0; i < blocksize; ++i) buffer1[i] += buffer2[i];

    out.put_block(current, blocksize, buffer1); 
  }

  delete[] buffer1;
  delete[] buffer2;
  return out;
};


template<class T>
PMODiagFile<T> PMODiagFile<T>::operator-(const PMODiagFile<T>& other) const {
  const int k = this->K_;
  const size_t fsize = this->filesize_;
  const int kk4 = std::max(k * k * 4, 1);
  const size_t blocksize = fsize / kk4;
  T* buffer1 = new T[blocksize];
  T* buffer2 = new T[blocksize];
  assert(fsize % kk4 == 0);

  PMODiagFile out(fsize, k, istart_, ifence_, jstart_, jfence_, false);

  size_t current = 0lu;
  for (int iblock = 0; iblock != kk4; ++iblock, current += blocksize) {
    this->get_block(current, blocksize, buffer1); 
    other.get_block(current, blocksize, buffer2);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize; ++i) buffer1[i] -= buffer2[i];

    out.put_block(current, blocksize, buffer1); 
  }

  delete[] buffer1;
  delete[] buffer2;
  return out;
};


template<class T>
T PMODiagFile<T>::get_energy_one_amp() const {
  const int isize = ifence_ - istart_;  
  const int jsize = jfence_ - jstart_;  
  const size_t ijsize = isize * jsize;
  assert(isize == jsize);
  const int k = this->K_;
  T* buffer = new T[ijsize * ijsize];
  T* coeff  = new T[ijsize * ijsize];

  T* c_coeff = coeff;
  // inefficient code
  const double geminal_exp = GEMINAL_EXP;
  for (int i1 = 0; i1 != isize; ++i1) {
    for (int i2 = 0; i2 != isize; ++i2) {
      for (int j1 = 0; j1 != isize; ++j1) {
        for (int j2 = 0; j2 != isize; ++j2, ++c_coeff) {
          if (i1 == j1 && i2 == j2 && i1 == i2) *c_coeff = 0.5 / (-geminal_exp);
          else if (i1 == j1 && i2 == j2) *c_coeff = 0.375 / (-geminal_exp);
          else if (i1 == j2 && i2 == j1) *c_coeff = 0.125 / (-geminal_exp);
          else *c_coeff = 0.0;
        }
      }
    }
  }

  int iblock = 0;
  T en = 0.0;
  for (int ki = -k; ki != std::max(k, 1); ++ki) {
    for (int kj = -k; kj != std::max(k, 1); ++kj, ++iblock) {
      this->get_block(iblock * ijsize * ijsize, ijsize * ijsize, buffer); 
//    #pragma omp parallel for reduction(+: en) schedule(dynamic, 100)
      for (int i1 = 0; i1 != isize; ++i1) {
        for (int i2 = 0; i2 != isize; ++i2) {
          const int xi1 = i2 + i1 * isize;
          const int xi2 = i1 + i2 * isize;
          for (int j1 = 0; j1 != isize; ++j1) {
            for (int j2 = 0; j2 != isize; ++j2) {
              const int xj = j2 + j1 * isize;
              en += buffer[xi1 * ijsize + xj] * (coeff[xi1 * ijsize + xj] * 2.0 - coeff[xi2 * ijsize + xj]);
            }
          }
        }
      }
    }
  }

  delete[] buffer;
  delete[] coeff;

  return en;
};


template<class T>
void PMODiagFile<T>::print() const {
  const int isize = ifence_ - istart_;  
  const int jsize = jfence_ - jstart_;  
  const size_t ijsize = isize * jsize;
  const int k = this->K_;
  T* buffer = new T[ijsize * ijsize];

  // I will print out in ne case...
  if (isize * jsize <= 25) {
    int iall = 0;
    for (int ki = -k; ki != std::max(k, 1); ++ki) {
      for (int kj = -k; kj != std::max(k, 1); ++kj, ++iall) {
        std::cout << "  block " << ki << " : " << kj << std::endl;
        this->get_block(iall * ijsize * ijsize, ijsize * ijsize, buffer); 
        const T* cbuf = buffer;
        for (int i = 0; i != ijsize; ++i) {
          for (int j = 0; j != ijsize; ++j, ++cbuf) {
            std::cout << std::setprecision(3) << *cbuf;
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
  delete[] buffer;
};


template<class T>
void PMODiagFile<T>::rprint() const {
  const int isize = ifence_ - istart_;  
  const int jsize = jfence_ - jstart_;  
  const size_t ijsize = isize * jsize;
  const int k = this->K_;
  T* buffer = new T[ijsize * ijsize];

  // I will print out in Ne or H2O cases...
  // singlet
  if (isize * jsize <= 25) {
    int iall = 0;
    for (int ki = -k; ki != std::max(k, 1); ++ki) {
      for (int kj = -k; kj != std::max(k, 1); ++kj, ++iall) {
        std::cout << "  block " << ki << " : " << kj << std::endl;
        this->get_block(iall * ijsize * ijsize, ijsize * ijsize, buffer); 
        const T* cbuf = buffer;
        for (int i = 0; i != ijsize; ++i) {
          for (int j = 0; j != ijsize; ++j, ++cbuf) {
            std::cout << std::setprecision(5) << std::setw(9) << (*cbuf).real();
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
  // triplet
  if (isize * jsize <= 25) {
    int iall = 0;
    for (int ki = -k; ki != std::max(k, 1); ++ki) {
      for (int kj = -k; kj != std::max(k, 1); ++kj, ++iall) {
        std::cout << "  block " << ki << " : " << kj << std::endl;
        this->get_block(iall * ijsize * ijsize, ijsize * ijsize, buffer); 
        for (int i1 = 0; i1 < isize; ++i1) {
          for (int i2 = 0; i2 < i1; ++i2) {
            const int xi1 = i2 + i1 * isize;
            const int xi2 = i1 + i2 * isize;
            for (int j1 = 0; j1 < jsize; ++j1) {
              for (int j2 = 0; j2 < j1; ++j2) {
                 const int xj = j2 + j1 * jsize;
                std::cout << std::setprecision(5) << std::setw(9) << (buffer[xi1 * ijsize + xj] - buffer[xi2 * ijsize + xj]).real();
              }
            }
            std::cout << std::endl;
          }
        }
        std::cout << std::endl;
      }
    }
  }
  delete[] buffer;
};


#endif
