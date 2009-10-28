//
// Author : Toru Shiozaki
// Date   : September 2009
//

#pragma once
#ifndef __src_util_pmofile_h
#define __src_util_pmofile_h

#include <src/pscf/pgeometry.h>
#include <src/util/pfile.h>
#include <src/pscf/f77.h>
#include <map>

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

    const boost::shared_ptr<PGeometry> geom_;
    size_t blocksize_;
    std::map<size_t, size_t> offset_;

  public:
    PMOFile(const boost::shared_ptr<PGeometry>,
        const long, const int, const int, const int,
        const int, const int,
        const int, const int,
        const int, const int,
        const bool late_init = false);
    PMOFile(const boost::shared_ptr<PGeometry>,
        const long, const int, const int, const int,
        const int, const int,
        const int, const int,
        const int, const int,
        const std::map<size_t, size_t>&, const bool late_init = false);
    ~PMOFile(); 

    std::map<size_t, size_t>::const_iterator oiter(int i, int j, int a, int b) const {
      const int k = this->K_;
      const int kk = k * 2;
      assert((i+j-a-b) % std::max(kk,1) == 0);
      const size_t tag = a+k+kk*(j+k+kk*(b+k));
      std::map<size_t, size_t>::const_iterator iter = offset_.find(tag);
      assert(iter != offset_.end());
      return iter;
    };
    void get_block2(int i, int j, int a, int b, T* data) const {
      this->get_block(oiter(i,j,a,b)->second, blocksize_, data);
    };
    void put_block2(int i, int j, int a, int b, const T* data) {
      this->put_block(oiter(i,j,a,b)->second, blocksize_, data);
    };
    void add_block2(int i, int j, int a, int b, const T* data) {
      this->add_block(oiter(i,j,a,b)->second, blocksize_, data);
    };

    void sort_inside_blocks();
    PMOFile<T> operator+(const PMOFile<T>&) const;
    PMOFile<T>& operator+=(const PMOFile<T>&);
    PMOFile<T> operator-(const PMOFile<T>&) const;
    PMOFile<T> operator*(const std::complex<double>&) const;

    boost::shared_ptr<PMOFile<T> > copy() const;
    void scale(double a);

    // Flip and add to this.
    void flip_symmetry();

    boost::shared_ptr<PMOFile<T> > contract(boost::shared_ptr<PMOFile<T> >, std::string);

    boost::shared_ptr<PMatrix1e> contract_density_J() const;
    boost::shared_ptr<PMatrix1e> contract_density_K() const;

    void print() const;
    void rprint() const;
    T get_energy_one_amp() const;

};


template<class T>
PMOFile<T>::PMOFile(const boost::shared_ptr<PGeometry> gm,
                    const long fsize, const int k,
                    const int istrt, const int ifen,
                    const int jstrt, const int jfen,
                    const int astrt, const int afen,
                    const int bstrt, const int bfen,
                    const bool late_init)
 : PFile<T>(fsize, k, late_init), geom_(gm),
   istart_(istrt), ifence_(ifen), jstart_(jstrt), jfence_(jfen),
   astart_(astrt), afence_(afen), bstart_(bstrt), bfence_(bfen) {

// current convention...
  blocksize_ = fsize / std::max(k*k*k*8, 1);
  size_t tag = 0;
  size_t current = 0;
  for (int b = -k; b != std::max(k,1); ++b) {
    for (int j = -k; j != std::max(k,1); ++j) {
      for (int a = -k; a != std::max(k,1); ++a, ++tag, current += blocksize_) {
        offset_.insert(std::make_pair(tag, current));
      }
    }
  }

};


template<class T>
PMOFile<T>::PMOFile(const boost::shared_ptr<PGeometry> gm,
                    const long fsize, const int k,
                    const int istrt, const int ifen,
                    const int jstrt, const int jfen,
                    const int astrt, const int afen,
                    const int bstrt, const int bfen,
                    const std::map<size_t, size_t>& of, const bool late_init)
 : PFile<T>(fsize, k, late_init), geom_(gm),
   istart_(istrt), ifence_(ifen), jstart_(jstrt), jfence_(jfen),
   astart_(astrt), afence_(afen), bstart_(bstrt), bfence_(bfen), offset_(of) {

// current convention...
  blocksize_ = fsize / std::max(k*k*k*8, 1);

};


template<class T>
PMOFile<T>::~PMOFile() {
};


// XXX conflicts in notation
// sort to chemist notations... but outermost block loops remain unchanged...
// that is unfavorable and needs some fix in PCompFile::mo_transform...
template<class T>
void PMOFile<T>::sort_inside_blocks() {

  T* buffer1 = new T[blocksize_];
  T* buffer2 = new T[blocksize_];

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const size_t absize = asize * bsize;

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);
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
    this->put_block(iter->second, blocksize_, buffer2);
  }

  delete[] buffer1;
  delete[] buffer2;
};


template<class T>
boost::shared_ptr<PMOFile<T> > PMOFile<T>::copy() const {
  boost::shared_ptr<PMOFile<T> > out(new PMOFile<T>(geom_, this->filesize_, this->K_,
                                                    istart_, ifence_, jstart_, jfence_,
                                                    astart_, afence_, bstart_, bfence_,
                                                    offset_, false));

  T* buffer = new T[blocksize_];
  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer);
    out->put_block(iter->second, blocksize_, buffer);
  }
  delete[] buffer;
  return out;
}


template<class T>
void PMOFile<T>::flip_symmetry() {
  const int k = this->K_;

  boost::shared_ptr<PMOFile<T> > other = this->copy();

  T* buffer1 = new T[blocksize_];
  T* buffer2 = new T[other->blocksize_];

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;

  assert(asize == bsize && isize == jsize);

  for (int kb = -k; kb != std::max(k, 1); ++kb) {
    for (int kj = -k; kj != std::max(k, 1); ++kj) {
      for (int ka = -k; ka != std::max(k, 1); ++ka) {
        int ki = ka+kb-kj;
        if (ki < -k) ki += k*2;
        else if (ki >=  k) ki -= k*2;

        get_block2(ki, kj, ka, kb, buffer1);
        other->get_block2(kj, ki, kb, ka, buffer2);

        T* cbuf = buffer1;
        for (int i = 0; i < isize; ++i) {
          for (int j = 0; j < jsize; ++j) {
            for (int a = 0; a < asize; ++a) {
              for (int b = 0; b < bsize; ++b, ++cbuf) {
                *cbuf += buffer2[a+asize*(b+bsize*(i+isize*j))];
              }
            }
          }
        }
        put_block2(ki, kj, ka, kb, buffer1);
      }
    }
  }

  delete[] buffer1;
  delete[] buffer2;
};


template<class T>
boost::shared_ptr<PMOFile<T> > PMOFile<T>::contract(boost::shared_ptr<PMOFile<T> > other, std::string jobname) {

  std::cout << "  Entering " << jobname << " contraction..." << std::endl;
  const int k = this->K_;
  const int kkk8 = std::max(k*k*k*8, 1);
  const size_t blocksize1 = blocksize_;
  const size_t blocksize2 = other->blocksize_;

  T* buffer1 = new T[blocksize1];
  T* buffer2 = new T[blocksize2];

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const int mfence = other->ifence_;
  const int mstart = other->istart_;
  const int msize = mfence - mstart;
  const int nfence = other->jfence_;
  const int nstart = other->jstart_;
  const int nsize = nfence - nstart;
  const int ijsize = isize * jsize;
  const int absize = asize * bsize;
  const int mnsize = msize * nsize;

  assert(blocksize1 == ijsize * absize && blocksize2 == mnsize * absize);

  T* target = new T[ijsize * mnsize];

  const size_t out_blocksize = ijsize * mnsize;
  const size_t out_filesize = kkk8 * out_blocksize;
  boost::shared_ptr<PMOFile<T> > out(new PMOFile<T>(geom_, out_filesize, k, istart_, ifence_, jstart_, jfence_,
                                                    mstart, mfence, nstart, nfence, false));

  const int k2 = std::max(k+k, 1);

  size_t current = 0lu; 
  for (int kj = -k; kj < std::max(k, 1); ++kj) { 
    for (int ki = -k; ki < std::max(k, 1); ++ki) {
      for (int kn = -k; kn < std::max(k, 1); ++kn) { 
        for (int km = -k; km < std::max(k, 1); ++km) {
          if ((km + kn - ki - kj) % k2 != 0) continue; 

          // block index is in Physicists' notation!!!
          const T zero = 0.0;
          fill(target, target + ijsize * mnsize, zero);

          for (int ka = -k; ka < std::max(k, 1); ++ka) {
            for (int kb = -k; kb < std::max(k, 1); ++kb) {
              if ((ka + kb - ki - kj) % k2 != 0) continue; 

              get_block2(kj, ki, kb, ka, buffer1);
              other->get_block2(kn, km, kb, ka, buffer2);
              const T one = 1.0;
              const T prefac = one / static_cast<T>(k2);
              zgemm_("C", "N", &ijsize, &mnsize, &absize, &prefac, buffer1, &absize, buffer2, &absize, &one, target, &ijsize);
            }
          }
          out->put_block2(kn, km, kj, ki, target);
        }
      }
    }
  }

  delete[] buffer1;
  delete[] buffer2;

  delete[] target;

  std::cout << "  done" << std::endl << std::endl;

  return out;
};


template<class T>
PMOFile<T> PMOFile<T>::operator-(const PMOFile<T>& other) const {

  assert(offset_.size() == other.offset_.size());
  T* buffer1 = new T[blocksize_];
  T* buffer2 = new T[blocksize_];

  PMOFile out(geom_, this->filesize_, this->K_, istart_, ifence_, jstart_, jfence_,
                                                astart_, afence_, bstart_, bfence_, offset_, false);

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);
    other.get_block(iter->second, blocksize_, buffer2);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize_; ++i) buffer1[i] -= buffer2[i];

    out.put_block(iter->second, blocksize_, buffer1);
  }

  delete[] buffer1;
  delete[] buffer2;
  return out;
};


template<class T>
PMOFile<T> PMOFile<T>::operator+(const PMOFile<T>& other) const {

  assert(offset_.size() == other.offset_.size());
  T* buffer1 = new T[blocksize_];
  T* buffer2 = new T[blocksize_];

  PMOFile out(geom_, this->filesize_, this->K_, istart_, ifence_, jstart_, jfence_,
                                                astart_, afence_, bstart_, bfence_, offset_, false);

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);
    other.get_block(iter->second, blocksize_, buffer2);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize_; ++i) buffer1[i] += buffer2[i];

    out.put_block(iter->second, blocksize_, buffer1);
  }

  delete[] buffer1;
  delete[] buffer2;
  return out;
};


template<class T>
PMOFile<T>& PMOFile<T>::operator+=(const PMOFile<T>& other) {

  assert(offset_.size() == other.offset_.size());
  T* buffer1 = new T[blocksize_];
  T* buffer2 = new T[blocksize_];

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);
    other.get_block(iter->second, blocksize_, buffer2);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize_; ++i) buffer1[i] += buffer2[i];

    this->put_block(iter->second, blocksize_, buffer1);
  }

  delete[] buffer1;
  delete[] buffer2;
  return *this;
};


template<class T>
PMOFile<T> PMOFile<T>::operator*(const std::complex<double>& a) const {

  T* buffer1 = new T[blocksize_];
  PMOFile out(geom_, this->filesize_, this->K_, istart_, ifence_, jstart_, jfence_,
                                                astart_, afence_, bstart_, bfence_, offset_, false);

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize_; ++i) buffer1[i] *= a;

    out.put_block(iter->second, blocksize_, buffer1);
  }

  delete[] buffer1;
  return out;
};


template<class T>
void PMOFile<T>::scale(double a) {

  T* buffer1 = new T[blocksize_];

  for (std::map<size_t, size_t>::const_iterator iter = offset_.begin(); iter != offset_.end(); ++iter) {
    this->get_block(iter->second, blocksize_, buffer1);

    #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < blocksize_; ++i) buffer1[i] *= a;

    this->put_block(iter->second, blocksize_, buffer1);
  }

  delete[] buffer1;
};


template<class T>
void PMOFile<T>::print() const {

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const size_t ijsize = isize * jsize;
  const size_t absize = asize * bsize;
  const int k = this->K_;
  const int kk = std::max(k + k, 1);
  T* buffer = new T[ijsize * absize];

  // I will print out in ne case...
  if (isize * jsize <= 25) {
    int iall = 0;
    for (int ki = -k; ki != std::max(k, 1); ++ki) {
      for (int kj = -k; kj != std::max(k, 1); ++kj) {
        for (int ka = -k; ka != std::max(k, 1); ++ka, ++iall) {
          std::cout << std::endl;
          this->get_block(iall * ijsize * absize, ijsize * absize, buffer);
          const T* cbuf = buffer;
          for (int i = 0; i != ijsize; ++i) {
            for (int j = 0; j != absize; ++j, ++cbuf) {
              std::cout << std::setprecision(4) << std::setw(15) << *cbuf;
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }
    }
  }
  delete[] buffer;
};


template<class T>
void PMOFile<T>::rprint() const {
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
      for (int kj = -k; kj != std::max(k, 1); ++kj) {
        for (int ka = -k; ka != std::max(k, 1); ++ka, ++iall) {
          std::cout << "  block " << ki << " : " << kj << " : " << ka << std::endl;
          int kb = ki + kj - ka;
          if (kb < -k) kb += 2*k;
          else if (kb >= k) kb -= 2*k;
          get_block2(ki, kj, ka, kb, buffer);
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
  }
  delete[] buffer;
};


template<class T>
T PMOFile<T>::get_energy_one_amp() const {
  const int isize = ifence_ - istart_;  
  const int jsize = jfence_ - jstart_;  
  const size_t ijsize = isize * jsize;
  assert(isize == jsize && isize == (afence_ - astart_) && jsize == (bfence_ - astart_));
  const int k = this->K_;
  T* buffer = new T[ijsize * ijsize];
  assert(geom_->gamma() > 0.0);
  const double gamma = geom_->gamma();
  const double ciiii = - 0.5 / gamma;
  const double cijij = - 0.375 / gamma;
  const double cijji = - 0.125 / gamma;
  int iblock = 0;
  T en = 0.0;
  for (int kb = -k; kb != std::max(k, 1); ++kb) {
    for (int kj = -k; kj != std::max(k, 1); ++kj) {
      for (int ka = -k; ka != std::max(k, 1); ++ka, ++iblock) {
        if (kb == kj && kb == ka) {
          this->get_block(iblock * ijsize * ijsize, ijsize * ijsize, buffer);
          for (int i1 = 0; i1 != isize; ++i1) { // i
            for (int i2 = 0; i2 != isize; ++i2) { // j
              if (i1 == i2) {
                const int xi = i2 + i1 * isize;
                const int xj = i2 + i1 * isize;
                en += buffer[xi * ijsize + xj] * ciiii;
              } else {
                const int xi = i2 + i1 * isize;
                const int xj1 = i2 + i1 * isize;
                const int xj2 = i1 + i2 * isize;
                en += buffer[xi * ijsize + xj1] * (2 * cijij - cijji);
                en += buffer[xi * ijsize + xj2] * (2 * cijji - cijij);
              }
            }
          }
        } else if (kb == kj && kb != ka) {
          this->get_block(iblock * ijsize * ijsize, ijsize * ijsize, buffer);
          for (int i1 = 0; i1 != isize; ++i1) { // i
            for (int i2 = 0; i2 != isize; ++i2) { // j
              const int xi = i2 + i1 * isize;
              const int xj = i2 + i1 * isize;
              en += buffer[xi * ijsize + xj] * (2 * cijij - cijji);
            }
          }
        } else if (ka == kj && kb != ka) {
          this->get_block(iblock * ijsize * ijsize, ijsize * ijsize, buffer);
          for (int i1 = 0; i1 != isize; ++i1) { // i
            for (int i2 = 0; i2 != isize; ++i2) { // j
              const int xi = i2 + i1 * isize;
              const int xj = i1 + i2 * isize;
              en += buffer[xi * ijsize + xj] * (2 * cijji - cijij);
            }
          }
        }
      }
    }
  }
  delete[] buffer;
  en /= static_cast<T>(std::max(k * k * 4, 1));
  return en;
};


// PMatrix1e assumes complex<double>.
template <class T>
boost::shared_ptr<PMatrix1e> PMOFile<T>::contract_density_J() const {
  typedef boost::shared_ptr<PMatrix1e> RefMatrix;

  // designed for operations like:
  // d^j_i * v^ip_jq -> h^p_q
  // d is assumed to be diagonal in Bloch orbitals.

  const int nocc = geom_->nocc() / 2;

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const size_t ijsize = isize * jsize;
  const size_t absize = asize * bsize;

  RefMatrix out(new PMatrix1e(geom_, bsize, jsize));

  const int k = this->K_;
  std::complex<double>* buffer = new std::complex<double>[blocksize_];

  for (int kb = -k; kb != std::max(k, 1); ++kb) {
    for (int kj = -k; kj != std::max(k, 1); ++kj) {
      for (int ka = -k; ka != std::max(k, 1); ++ka) {

        int ki = ka+kb-kj;
        if (ki < -k) ki += k*2;
        else if (ki >= k) ki -= k*2;

        if (ki != ka) continue;
        if (kb != kj) continue;

        std::complex<double>* oblock = out->bpw(kb);

        get_block2(ki, kj, ka, kb, buffer);
        const std::complex<double>* cbuf = buffer;
        // Assumes that it is already sorted.
        for (int i = 0; i != isize; ++i) {
          for (int j = 0; j != jsize; ++j) {
            for (int a = 0; a != asize; ++a) {
              const double cden = (i == a && i+istart_ <= nocc) ? 2.0/std::max(k+k,1) : 0.0;
              for (int b = 0; b != bsize; ++b, ++cbuf) {
                oblock[b + bsize*j] += *cbuf * cden;
              }
            }
          }
        }
      }
    }
  }

  delete[] buffer;
  return out;
};


// PMatrix1e assumes complex<double>.
template <class T>
boost::shared_ptr<PMatrix1e> PMOFile<T>::contract_density_K() const {
  typedef boost::shared_ptr<PMatrix1e> RefMatrix;

  // designed for operations like:
  // d^j_i * v^ip_qj -> K^p_q
  // d is assumed to be diagonal in Bloch orbitals.

  const int nocc = geom_->nocc() / 2;

  const int isize = ifence_ - istart_;
  const int jsize = jfence_ - jstart_;
  const int asize = afence_ - astart_;
  const int bsize = bfence_ - bstart_;
  const size_t ijsize = isize * jsize;
  const size_t absize = asize * bsize;

  RefMatrix out(new PMatrix1e(geom_, asize, jsize));

  const int k = this->K_;
  std::complex<double>* buffer = new std::complex<double>[blocksize_];

  for (int kb = -k; kb != std::max(k, 1); ++kb) {
    for (int kj = -k; kj != std::max(k, 1); ++kj) {
      for (int ka = -k; ka != std::max(k, 1); ++ka) {

        int ki = ka+kb-kj;
        if (ki < -k) ki += k*2;
        else if (ki >= k) ki -= k*2;

        if (ki != kb) continue;
        if (ka != kj) continue;

        std::complex<double>* oblock = out->bpw(ka);

        get_block2(ki, kj, ka, kb, buffer);
        const std::complex<double>* cbuf = buffer;
        // Assumes that it is already sorted.
        for (int i = 0; i != isize; ++i) {
          for (int j = 0; j != jsize; ++j) {
            for (int a = 0; a != asize; ++a) {
              for (int b = 0; b != bsize; ++b, ++cbuf) {
                const double cden = (i == b && i+istart_ <= nocc) ? 1.0/std::max(k+k,1) : 0.0;
                oblock[a + asize*j] += *cbuf * cden;
              }
            }
          }
        }
      }
    }
  }

  delete[] buffer;
  return out;
};


#endif
