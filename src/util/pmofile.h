//
// Author : Toru Shiozaki
// Date   : September 2009
//

#pragma once
#ifndef __src_util_pmofile_h
#define __src_util_pmofile_h

#include <src/util/pfile.h>
#include <src/pscf/f77.h>
#include <src/util/pmodiagfile.h>

template<class T>
class PMOFile : public PFile<T> {
  protected:
    const int istart_;
    const int jstart_;
    const int astart_;
    const int bstart_;
    const int ifence_;
    const int jfence_;
    const int afence_;
    const int bfence_;

  public:
    PMOFile(const long, const int, const int, const int,
                                   const int, const int,
                                   const int, const int,
                                   const int, const int, const bool late_init = false); 
    ~PMOFile(); 

    void sort_inside_blocks();

    boost::shared_ptr<PMODiagFile<T> > contract(boost::shared_ptr<PMOFile<T> >, std::string);
};


template<class T>
PMOFile<T>::PMOFile(const long fsize, const int k,
                    const int istrt, const int ifen,
                    const int jstrt, const int jfen,
                    const int astrt, const int afen,
                    const int bstrt, const int bfen, const bool late_init)
 : PFile<T>(fsize, k, late_init), istart_(istrt), ifence_(ifen), jstart_(jstrt), jfence_(jfen),
                                  astart_(astrt), afence_(afen), bstart_(bstrt), bfence_(bfen) {


};


template<class T>
PMOFile<T>::~PMOFile() {
};


// sort to chemist notations... but outermost block loops remain unchanged...
// that is unfavorable and needs some fix in PCompFile::mo_transform... TODO
template<class T>
void PMOFile<T>::sort_inside_blocks() {
  const int k = this->K_;
  const int KKK8 = k * k * k * 8;
  assert(this->filesize_ % KKK8 == 0);

  const size_t blocksize = this->filesize_ / KKK8;
  size_t current = 0lu; 

  T* buffer1 = new T[blocksize];
  T* buffer2 = new T[blocksize];

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const size_t absize = asize * bsize;

  for (int iblock = 0; iblock != KKK8; ++iblock, current += blocksize) { 
    this->get_block(current, blocksize, buffer1); 
    #pragma omp parallel for
    for (int i = 0; i < isize; ++i) {
      const T* cbuf = buffer1 + bsize * asize * jsize * i;
      for (int a = 0; a != asize; ++a) {
        T* cbuf2 = buffer2 + bsize * (a + asize * (0 + jsize * i));
        for (int j = 0; j != jsize; ++j, cbuf += bsize, cbuf2 += absize) {
          ::memcpy(cbuf2, cbuf, bsize * sizeof(T));
        }
      }
    }
    this->put_block(current, blocksize, buffer2);
  }

  delete[] buffer1;
  delete[] buffer2;
};


template<class T>
boost::shared_ptr<PMODiagFile<T> > PMOFile<T>::contract(boost::shared_ptr<PMOFile<T> > other, std::string jobname) {

  std::cout << "  Entering " << jobname << " contraction..." << std::endl;
  const int k = this->K_;
  const int KKK8 = std::max(k * k * k * 8, 1);
  assert(this->filesize_ % KKK8 == 0);
  const size_t blocksize = this->filesize_ / KKK8;

  T* buffer1 = new T[blocksize];
  T* buffer2 = new T[blocksize];

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const int ijsize = isize * jsize;
  const int absize = asize * bsize;
  assert(ijsize * absize == blocksize);

  T* target = new T[ijsize * ijsize];

  const size_t out_blocksize = ijsize * ijsize;
  const size_t out_filesize = std::max(k * k * 4, 1) * out_blocksize; 
  boost::shared_ptr<PMODiagFile<T> > out(new PMODiagFile<T>(out_filesize, k, istart_, ifence_, jstart_, jfence_, false));

  size_t current = 0lu; 
  for (int kb = -k; kb < std::max(k, 1); ++kb) { 
    for (int kj = -k; kj < std::max(k, 1); ++kj) { 
      for (int ka = -k; ka < std::max(k, 1); ++ka, current += blocksize) { 
        int ki = ka + kb - kj; 
        if (ki < - k) ki += k * 2; 
        else if (ki >=  k) ki -= k * 2; 

        this->get_block(current, blocksize, buffer1);
        other->get_block(current, blocksize, buffer2);

        const T one = 1.0;
        const T zero = 0.0;

        // TODO how to deal with this???
        // assuming complex<double>....
        zgemm_("C", "N", &ijsize, &ijsize, &absize, &one, buffer1, &absize, buffer2, &absize, &zero, target, &ijsize);

        const size_t kij = (ki + k) + (kj + k) * std::max(k + k, 1);
        out->add_block(kij * out_blocksize, out_blocksize, target);

      }
    }
  }

  delete[] buffer1;
  delete[] buffer2;

  delete[] target;

  std::cout << "  done" << std::endl << std::endl;

  return out;
};


#endif

