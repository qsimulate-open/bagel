//
// Author : Toru Shiozaki
// Date   : September 2009
//

#ifndef __src_util_pmofile_h
#define __src_util_pmofile_h

#include <src/util/pfile.h>

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
  const size_t bsize = bfence_ - bstart_;
  const size_t absize = asize * bsize;

  for (int iblock = 0; iblock != KKK8; ++iblock, current += blocksize) { 
    get_block(current, blocksize, buffer1); 
    const T* cbuf = buffer1;
    #pragma omp parallel for
    for (int i = 0; i < isize; ++i) {
      for (int a = 0; a != asize; ++a) {
        T* cbuf2 = buffer2 + bsize * (a + asize * jsize * i);
        for (int j = 0; j != jsize; ++j, cbuf += bsize, cbuf2 += absize)
          ::memcpy(cbuf2, cbuf, bsize * sizeof(T));
      }
    }
    put_block(current, blocksize, buffer2);
  }

  delete[] buffer1;
  delete[] buffer2;
};


#endif

