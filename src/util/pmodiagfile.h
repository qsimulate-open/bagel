//
// Author : Toru Shiozaki
// Date   : September 2009
//

#pragma once
#ifndef __src_util_pmodiagfile_h
#define __src_util_pmodiagfile_h

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

    #pragma omp parallel for schedule(dynamic, 128)
    for (int i = 0; i < blocksize; ++i) buffer1[i] -= buffer2[i];

    out.put_block(current, blocksize, buffer1); 
  }

  delete[] buffer1;
  delete[] buffer2;
  return out;
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


#endif
