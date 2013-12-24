//
// BAGEL - Parallel electron correlation program.
// Filename: pcompfile.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// compressed file of double

#ifndef __src_util_pcompfile_h
#define __src_util_pcompfile_h

#include <stddef.h>
#include <cstring>
#include <cassert>
#include <vector>
#include <fstream>
#include <memory>
#include <iostream>
#include <algorithm>
#include <src/util/constants.h>
#include <src/pscf/pgeometry.h>
#include <src/pscf/pcoeff.h>
#include <src/pscf/pmofile.h>
#include <src/pscf/filename.h>
#include <src/pscf/cache.h>
#include <src/util/f77.h>
#include <src/math/algo.h>

namespace bagel {

template<class T>
class PCompFile {
  protected:
    const std::shared_ptr<PGeometry> geom_;

    // Exponent of STG for this file. Can be different from geom_->gamma()
    const double gamma_;

    std::shared_ptr<std::fstream> file_;
    long filesize_;
    std::string filename_;
    std::string jobname_;
    std::vector<double> schwarz_;

    std::vector<size_t> num_int_each_;

    std::vector<std::shared_ptr<const Shell>> basis_;
    std::vector<int> offset_;
    size_t max_num_int_;

    int K_, L_, S_;
    double A_;

  public:
    PCompFile(std::shared_ptr<PGeometry>,  const double gam,
        const bool late_init = false, const std::string jobname = "");
    ~PCompFile();

    // create file...
    const std::string filename() const { return filename_; };
    void append(const long, const double*);
    // evaluates one block of integrals.
    virtual void eval_new_block(double*, int, int, int);

    // append mode cannot be read, hence reopen. Add, Put and Get operation has been implemented
    void reopen_with_inout();
    void add_block(const long, const long, const double*);
    void put_block(const long, const long, const double*);
    void get_block(const long, const long, double*);

    // returns some misc infos
    double schwarz(int i) const { return schwarz_[i]; };
    std::vector<int> offset() const { return offset_; };
    int offset(size_t i) const { return offset_[i]; };
    std::vector<std::shared_ptr<const Shell>> basis() const { return basis_; };
    std::shared_ptr<const Shell> basis(size_t i) const { return basis_[i]; };
    size_t nbasis(size_t i) const { return basis_[i]->nbasis(); };
    int basissize() const { return basis_.size(); };

    // calculate integrals here
    void init_schwarz();
    std::vector<size_t> num_int_each() const { return num_int_each_; };
    size_t num_int_each(const int i) const { return num_int_each_[i]; };
    // will be redefined in the derived classes
    virtual void calculate_num_int_each();
    virtual void store_integrals();

    // or import integrals computed externally
    void set_schwarz(const std::vector<double>& a) { schwarz_ = a; };
    void set_num_int_each(const std::vector<size_t>& nie) { num_int_each_ = nie;
                                                            max_num_int_ = *std::max_element(nie.begin(), nie.end()); };
    void set_integrals(std::shared_ptr<std::fstream> inp) { file_ = inp; };

    size_t max_num_int() const { assert(max_num_int_ > 0lu); return max_num_int_; };

    int K() const { return K_; };
    int L() const { return L_; };
    int S() const { return S_; };
    double A() const { return A_; };

    // AO-to-MO integral transformation...
    std::shared_ptr<PMOFile<std::complex<double>>>
                mo_transform(std::shared_ptr<PCoeff>,
                             std::shared_ptr<PCoeff>,
                             std::shared_ptr<PCoeff>,
                             std::shared_ptr<PCoeff>,
                             const int istart, const int ifence,
                             const int jstart, const int jfence,
                             const int astart, const int afence,
                             const int bstart, const int bfence,
                             const std::string jobname = "intermediate");

};


template<class T>
PCompFile<T>::PCompFile(std::shared_ptr<PGeometry> gm, const double gam, const bool late_init, const std::string jobname)
 : geom_(gm), gamma_(gam), jobname_(jobname) {

  Filename tmpf;
  filename_ = tmpf.filename_next();

  std::shared_ptr<std::fstream> tmp(new std::fstream(filename_.c_str(), std::ios::out | std::ios::binary));
  file_ = tmp;

  { // prepare offset and basis
    typedef std::shared_ptr<const Atom> RefAtom;
    typedef std::shared_ptr<const Shell> RefShell;

    const std::vector<RefAtom> atoms = geom_->atoms();
    int cnt = 0;
    for (std::vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
      const std::vector<RefShell> tmp = (*aiter)->shells();
      basis_.insert(basis_.end(), tmp.begin(), tmp.end());
      const std::vector<int> tmpoff = geom_->offset(cnt);
      offset_.insert(offset_.end(), tmpoff.begin(), tmpoff.end());
    }
  }

  S_ = geom_->S();
  L_ = geom_->L();
  K_ = geom_->K();
  A_ = geom_->A();

  // late_init is true when constructed from PairCompFile class
  if (!late_init) {
    init_schwarz();
    calculate_num_int_each();
  }
};


template<class T>
PCompFile<T>::~PCompFile() {
  unlink(filename_.c_str());
};


template<class T>
void PCompFile<T>::add_block(const long position, const long length, const double* data) {
  long remaining = length;
  long current = 0L;
  const double cone = 1.0;
  double* work = (double*) work_char;

  while (remaining > 0L) {
    const long readsize = std::min(remaining, cachesize) * sizeof(double);
    const long rsize = std::min(remaining, cachesize);

    file_->clear();
    file_->seekg((position + current) * sizeof(double));
    file_->read((char*)work_char, readsize);
    for (long i = 0; i != rsize; ++i) work[i] += data[current + i];

    file_->clear();
    file_->seekp((position + current) * sizeof(double));
    file_->write((const char*)work_char, readsize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PCompFile<T>::get_block(const long position, const long length, double* data) {
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long readsize = std::min(remaining, cachesize) * sizeof(double);
    file_->clear();
    file_->seekg((position + current) * sizeof(double));
    file_->read((char*)(data + current), readsize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PCompFile<T>::put_block(const long position, const long length, const double* data) {
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long writesize = std::min(remaining, cachesize) * sizeof(double);
    file_->clear();
    file_->seekp((position + current) * sizeof(double));
    file_->write((const char*)(data + current), writesize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PCompFile<T>::append(const long length, const double* data) {
  file_->clear();
  file_->seekp(0, std::ios::end);
  long remaining = length;
  long current = 0L;

  while (remaining > 0L) {
    const long writesize = std::min(remaining, cachesize) * sizeof(double);
    file_->write((const char*)(data + current), writesize);

    remaining -= cachesize;
    current += cachesize;
  }
};


template<class T>
void PCompFile<T>::reopen_with_inout() {
  file_->close();
  file_->open(filename_.c_str(), std::ios::in | std::ios::out | std::ios::binary);
};


template<class T>
void PCompFile<T>::calculate_num_int_each() {

  typedef std::shared_ptr<const Shell> RefShell;

  unsigned long data_written = 0ul;
  num_int_each_.resize((S_ + S_ + 1) * (S_ + S_ + 1) * (L_ + 1));
  const int size = basis_.size(); // number of shells

//#pragma omp parallel for reduction (+:data_written)
  for (int m1 = - S_; m1 <= S_; ++m1) {
    const double m1disp[3] = {0.0, 0.0, m1 * A_};
    size_t offset = (m1 + S_) * (L_ + 1) * (S_ * 2 + 1);
    for (int m2 = 0; m2 <= L_; ++m2) { // use bra-ket symmetry!!!
      const double m2disp[3] = {0.0, 0.0, m2 * A_};
      for (int m3 = m2 - S_; m3 <= m2 + S_; ++m3, ++offset) {
        const double m3disp[3] = {0.0, 0.0, m3 * A_};
        size_t thisblock = 0ul;
        for (int i0 = 0; i0 != size; ++i0) {
          const int b0offset = offset_[i0];
          const int b0size = basis_[i0]->nbasis();

          for (int i1 = 0; i1 != size; ++i1) {
            const int b1offset = offset_[i1];
            const int b1size = basis_[i1]->nbasis();

            for (int i2 = 0; i2 != size; ++i2) {
              const int b2offset = offset_[i2];
              const int b2size = basis_[i2]->nbasis();

              for (int i3 = 0; i3 != size; ++i3) {
                const int b3offset = offset_[i3];
                const int b3size = basis_[i3]->nbasis();

                const double integral_bound = schwarz_[(m1 + K_) * size * size + i0 * size + i1]
                                            * schwarz_[(m3 - m2 + K_) * size * size + i2 * size + i3];
                const bool skip_schwarz = integral_bound < schwarz_thresh__;
                if (skip_schwarz) continue;
                data_written += b0size * b1size * b2size * b3size;
                thisblock += b0size * b1size * b2size * b3size;
              }
            }
          }
        }
        num_int_each_[offset] = thisblock;
      }
    }
  }

  max_num_int_ = *std::max_element(num_int_each_.begin(), num_int_each_.end());

  if (jobname_ != "NULL") {
    std::cout << std::fixed << "  Using ";
    const size_t data_written_byte = data_written * sizeof(double);
    if (data_written_byte > 1.0e9) {
      std::cout << std::setprecision(1) << data_written_byte / 1.0e9 << " GB";
    } else if (data_written_byte > 1.0e6) {
      std::cout << std::setprecision(1) << data_written_byte / 1.0e6 << " MB";
    } else {
      std::cout << std::setprecision(1) << data_written_byte / 1.0e3 << " KB";
    }

    std::cout << " hard disk for storing \"" << jobname_ << "\"" << std::endl;
  }
};


template<class T>
void PCompFile<T>::eval_new_block(double* out, int m1, int m2, int m3) {

  typedef std::shared_ptr<const Shell> RefShell;

  const double m1disp[3] = {0.0, 0.0, m1 * A_};
  const double m2disp[3] = {0.0, 0.0, m2 * A_};
  const double m3disp[3] = {0.0, 0.0, m3 * A_};

  const int size = basis_.size(); // number of shells
  size_t* blocks = new size_t[size * size * size * size + 1];
  blocks[0] = 0;
  int iall = 0;
  for (int i0 = 0; i0 != size; ++i0) {
    const int b0size = basis_[i0]->nbasis();
    for (int i1 = 0; i1 != size; ++i1) {
      const int b1size = basis_[i1]->nbasis();
      for (int i2 = 0; i2 != size; ++i2) {
        const int b2size = basis_[i2]->nbasis();
        for (int i3 = 0; i3 != size; ++i3, ++iall) {
          const int b3size = basis_[i3]->nbasis();
          const double integral_bound = schwarz_[(m1 + K_) * size * size + i0 * size + i1]
                                      * schwarz_[(m3 - m2 + K_) * size * size + i2 * size + i3];
          const bool skip_schwarz = integral_bound < schwarz_thresh__;
          blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
        }
      }
    }
  }
//#pragma omp parallel for
  for (int i01 = 0; i01 < size * size; ++i01) {
    const int i1 = i01 % size;
    const int i0 = (i01 - i1) / size;
    {
      int offset = i01 * size * size;
      const RefShell b0 = basis_[i0]; // b0 is the center cell
      const RefShell b1 = basis_[i1]->move_atom(m1disp);
      for (int i2 = 0; i2 != size; ++i2) {
        const RefShell b2 = basis_[i2]->move_atom(m2disp);
        for (int i3 = 0; i3 < size; ++i3, ++offset) {
          const RefShell b3 = basis_[i3]->move_atom(m3disp);

          if (blocks[offset] == blocks[offset + 1]) continue;

          std::array<RefShell,4> input = {{ b3, b2, b1, b0 }};

          T batch(input, 1.0, gamma_);
          batch.compute();
          const double* bdata = batch.data();
          ::memcpy(out + blocks[offset], bdata, (blocks[offset + 1] - blocks[offset]) * sizeof(double));
        }
      }
    }
  }
  delete[] blocks;
};


template<class T>
void PCompFile<T>::store_integrals() {
  // TODO control the cache size
  const size_t cachesize_max = 120000000lu;
  const size_t cachesize = std::max(cachesize_max, max_num_int_);
  double* dcache = new double[cachesize];
  size_t remaining = cachesize;
  size_t current = 0lu;
  size_t cnt = 0lu;
  for (int m1 = -S_; m1 <= S_; ++m1) {
    for (int m2 = 0; m2 <= L_; ++m2) { // use bra-ket symmetry!!!
      for (int m3 = m2 - S_; m3 <= m2 + S_; ++m3, ++cnt) {
        if (remaining < num_int_each(cnt)) {
          append(current, dcache);
          current = 0lu;
          remaining = cachesize;
        }
        eval_new_block(&dcache[current], m1, m2, m3);
        current += num_int_each(cnt);
        remaining -= num_int_each(cnt);
      }
    }
  }
  append(current, dcache);
  delete[] dcache;
};


template<class T>
void PCompFile<T>::init_schwarz() {
  typedef std::shared_ptr<const Shell> RefShell;
  typedef std::shared_ptr<const Atom> RefAtom;

  const int size = basis_.size(); // the number of shells per unit cell
  schwarz_.resize(size * size * (2 * K_ + 1));

//#pragma omp parallel for
  for (int m = - K_; m <= K_; ++m) {
    const double disp[3] = {0.0, 0.0, m * A_};
    for (int i0 = 0; i0 != size; ++i0) { // center unit cell
      const RefShell b0 = basis_[i0];
      for (int i1 = 0; i1 != size; ++i1) {
        const RefShell b1 = basis_[i1]->move_atom(disp);

        std::array<RefShell,4> input = {{b0, b1, b0, b1}};
        T batch(input, 1.0, gamma_);
        batch.compute();
        const double* data = batch.data();
        const int datasize = batch.data_size();
        double cmax = 0.0;
        for (int xi = 0; xi != datasize; ++xi, ++data) {
          const double absed = std::fabs(*data);
          if (absed > cmax) cmax = absed;
        }
        const size_t boffset = static_cast<size_t>(m + K_) * size * size;
        schwarz_[boffset + i0 * size + i1] = std::sqrt(cmax);
      }
    }
  }
};


template<class T>
std::shared_ptr<PMOFile<std::complex<double>>>
  PCompFile<T>::mo_transform(std::shared_ptr<PCoeff> coeff_i,
                             std::shared_ptr<PCoeff> coeff_j,
                             std::shared_ptr<PCoeff> coeff_a,
                             std::shared_ptr<PCoeff> coeff_b,
                             const int istart, const int ifence,
                             const int jstart, const int jfence,
                             const int astart, const int afence,
                             const int bstart, const int bfence,
                             const std::string jobname) {
// Loading a (2K * 2K * nov) quantity on memory

  const int isize = ifence - istart;
  const int jsize = jfence - jstart;
  const int asize = afence - astart;
  const int bsize = bfence - bstart;
  const size_t noovv = static_cast<size_t>(isize) * jsize * asize * bsize;
  assert(noovv > 0);

  const double pi = 3.14159265358979323846;
  const int KK = K_ + K_;

  const int maxK1 = std::max(K_, 1);
  const std::complex<double> czero(0.0, 0.0);
  const std::complex<double> cone(1.0, 0.0);

  const int unit = 1;

  const int nbasis1 = geom_->nbasis();
  const int nbasis2 = nbasis1 * nbasis1;
  const int nbasis3 = nbasis2 * nbasis1;
  const size_t nbasis4 = static_cast<size_t>(nbasis2) * nbasis2;

  const size_t filesize = noovv * std::max(KK, 1) * std::max(KK, 1) * std::max(KK, 1);
  std::cout << "  Creating " << jobname << "  of size ";
  const size_t filesize_byte = filesize * sizeof(std::complex<double>);
  if (filesize_byte > 1.0e9) {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e9 << " GB" << std::endl;
  } else if (filesize_byte > 1.0e6) {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e6 << " MB" << std::endl;
  } else {
    std::cout << std::setprecision(1) << filesize_byte / 1.0e3 << " KB" << std::endl;
  }
  std::shared_ptr<PMOFile<std::complex<double>>>
    mo_int(new PMOFile<std::complex<double>>(geom_, filesize, K_,
                                              istart, ifence, jstart, jfence,
                                              astart, afence, bstart, bfence, true));

  // we are assuming that the (c.c. two-electron integrals for a unit cell)*K^2 can be
  // held in core. If that is not the case, this must be rewritten.

  const size_t nv = nbasis3 * bsize;
  const size_t nov = nbasis2 * jsize * bsize;
  const size_t novv = nbasis1 * jsize * asize * bsize;
  // allocating a temp array

  const size_t nmax = std::max(nbasis4, std::max((size_t)std::max(std::max(nv, nov), novv), noovv));
  std::complex<double>* data = new std::complex<double>[nmax];
  const size_t nalloc = std::max((size_t)nmax, std::max((size_t)std::max(nov,novv), noovv));

  std::complex<double>* datas = new std::complex<double>[nalloc];
  std::complex<double>* conjc = new std::complex<double>[nbasis1 * std::max(isize, jsize)];
  double* data_read = new double[nbasis4];

  const int size = basis_.size();
  int* blocks = new int[size * size * size * size + 1];

  const size_t sizem1 = S_ + S_ + 1lu;
  const size_t sizem2 = L_ + L_ + 1lu;
  const size_t sizem2abs = L_ + 1lu;
  const size_t num_loops = sizem1 * sizem1 * sizem2;
  size_t loop_counter = 0lu;
  size_t loop_mod10 = 0lu;

  size_t allocsize = *std::max_element(num_int_each_.begin(), num_int_each_.end());

  PFile<std::complex<double>> intermediate_mmK(nv * std::max(KK, 1), K_, false);
  std::complex<double>* intermediate_novv = new std::complex<double>[novv];
  PFile<std::complex<double>> intermediate_mKK(std::max(nov * KK * KK, nov), K_, false);
  PFile<std::complex<double>> intermediate_KKK(std::max(novv * KK * KK * KK, novv), K_, true);

  for (int q1 = -S_; q1 <= S_; ++q1) {
    const bool q1_front = (q1 == -S_);
    const int m1 = q1;
    intermediate_mKK.clear();
    for (int q2 = -L_; q2 <= L_; ++q2) {
      const int m2 = q2;
      intermediate_mmK.clear();
      for (int q3 = -S_; q3 <= S_; ++q3, ++loop_counter) {
        const int m3 = m2 + q3;

        const double* cdata;
        if (q2 >= 0) {
          const size_t key = q3 + S_ + sizem1 * (q2 + sizem2abs * (q1 + S_));
          size_t datasize_acc = 0lu;
          for (int i = 0; i != key; ++i) datasize_acc += num_int_each_[i];
          get_block(datasize_acc, num_int_each_[key], data_read);
          cdata = data_read;
        }
        else { // q2 < 0
          const size_t key = q1 + S_ + sizem1 * (- q2 + sizem2abs * (q3 + S_));
          size_t datasize_acc = 0lu;
          for (int i = 0; i != key; ++i) datasize_acc += num_int_each_[i];
          get_block(datasize_acc, num_int_each_[key], data_read);
          cdata = data_read;
        }

        if (loop_mod10 < loop_counter * 10lu / num_loops) {
          loop_mod10++;
          std::cout << "  Loop " << loop_mod10 * 10lu << " percent done" << std::endl;
        }

        size_t local_counter = 0lu;

        blocks[0] = 0;
        int iall = 0;
        for (int i0 = 0; i0 != size; ++i0) {
          const int b0size = basis_[i0]->nbasis();
          for (int i1 = 0; i1 != size; ++i1) {
            const int b1size = basis_[i1]->nbasis();
            for (int i2 = 0; i2 != size; ++i2) {
              const int b2size = basis_[i2]->nbasis();
              for (int i3 = 0; i3 != size; ++i3, ++iall) {
                const int b3size = basis_[i3]->nbasis();

                double integral_bound;
                if (m2 >= 0) {
                  integral_bound = schwarz_[((q1 + K_) * size + i0) * size + i1]
                                 * schwarz_[((q3 + K_) * size + i2) * size + i3];
                } else {
                  integral_bound = schwarz_[((q1 + K_) * size + i2) * size + i3]
                                 * schwarz_[((q3 + K_) * size + i0) * size + i1];
                }
                const bool skip_schwarz = integral_bound < schwarz_thresh__;
                blocks[iall + 1] = blocks[iall] + (skip_schwarz ? 0 : (b0size * b1size * b2size * b3size));
              }
            }
          }
        }

//      #pragma omp parallel for
        for (int i0 = 0; i0 < size; ++i0) {
          int noffset = i0 * size * size * size;
          const int b0offset = offset(i0);
          const int b0size = nbasis(i0);
          for (int i1 = 0; i1 != size; ++i1) {
            const int b1offset = offset(i1);
            const int b1size = nbasis(i1);
            for (int i2 = 0; i2 != size; ++i2) {
              const int b2offset = offset(i2);
              const int b2size = nbasis(i2);
              for (int i3 = 0; i3 != size; ++i3, ++noffset) {
                const int b3offset = offset(i3);
                const int b3size = nbasis(i3);

                if (blocks[noffset] != blocks[noffset + 1]) {
                  const double* ndata = cdata + blocks[noffset];
                  std::complex<double>* label = data + b3offset + nbasis1 * (b2offset + nbasis1 * (b1offset + nbasis1 * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += nbasis3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += nbasis2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += nbasis1) {
                        for (int j3 = 0; j3 != b3size; ++j3, ++label, ++ndata) {
                          *label = static_cast<std::complex<double>>(*ndata);
                        }
                        label -= b3size;
                      }
                      label -= b2size * nbasis1;
                    }
                    label -= b1size * nbasis2;
                  }
                } else {
                  std::complex<double>* label = data + b3offset + nbasis1 * (b2offset + nbasis1 * (b1offset + nbasis1 * b0offset));
                  for (int j0 = 0; j0 != b0size; ++j0, label += nbasis3) {
                    for (int j1 = 0; j1 != b1size; ++j1, label += nbasis2) {
                      for (int j2 = 0; j2 != b2size; ++j2, label += nbasis1) {
                        std::fill(label, label + b3size, czero);
                      }
                      label -= b2size * nbasis1;
                    }
                    label -= b1size * nbasis2;
                  }
                } // end of if skip_schwarz

              }
            }
          }
        } // end of shell loops; now data is ready

        if (m2 >= 0) { // start from j & b contractions -- needs transposition
          const int mn = nbasis2;
          blas::transpose(data, mn, mn, datas);
          ::memcpy(data, datas, nbasis4 * sizeof(std::complex<double>));
        }

//      #pragma omp parallel for
        for (int nkb = -K_; nkb < maxK1; ++nkb) {
          const int nkbc = nkb + K_;
          const std::complex<double> exponent(0.0, K_ != 0 ? (m3 * nkb * pi) / K_ : 0.0);
          const std::complex<double> prefac = exp(exponent);
          int offset1 = 0;
          int offset2 = 0;
          std::fill(datas, datas+nv, 0.0);
//        #pragma omp parallel for
          for (int ii = 0; ii < nbasis1; ++ii) {
            const size_t offset1 = nbasis3 * ii;
            const size_t offset2 = nbasis2 * bsize * ii;
            zgemm3m_("N", "N", &nbasis2, &bsize, &nbasis1, &prefac, data + offset1, &nbasis2,
                                                         coeff_b->bp(nkb) + nbasis1 * bstart, &nbasis1, &cone,
                                                         datas + offset2, &nbasis2);
          }
          intermediate_mmK.add_block(nv*nkbc, nv, datas);
        } // end of contraction b for given m3

      } // end of m3 loop

      // intermediate_mmK is ready
      for (int nkb = -K_; nkb < maxK1; ++nkb) {
        const int nkbc = nkb + K_;
        int nbj = nkbc * KK;
        intermediate_mmK.get_block(nv*nkbc, nv, data);
        for (int nkj = -K_; nkj < maxK1; ++nkj, ++nbj) {
          std::complex<double>* conjc2 = new std::complex<double>[nbasis1 * jsize];
          const std::complex<double> exponent(0.0, K_ != 0 ? (- m2 * nkj * pi) / K_ : 0.0);
          const std::complex<double> prefac = exp(exponent);

          const std::complex<double>* jdata = coeff_j->bp(nkj) + nbasis1 * jstart;
          for (int ii = 0; ii != nbasis1 * jsize; ++ii) conjc2[ii] = conj(jdata[ii]);
          const int nsize = nbasis2 * bsize;
          fill(datas, datas+nov, czero);
          zgemm3m_("N", "N", &nsize, &jsize, &nbasis1, &prefac, data, &nsize,
                                                              conjc2, &nbasis1, &cone,
                                                              datas, &nsize);
          delete[] conjc2;
          intermediate_mKK.add_block(nbj*nov, nov, datas);
        }
      } // end of contraction j for given m2

    } // end of m2 loop

    // intermediate_mKK is ready
    for (int nkb = -K_, nbja = 0, nbj = 0; nkb != maxK1; ++nkb) {
      for (int nkj = -K_; nkj != maxK1; ++nkj, ++nbj) {
        intermediate_mKK.get_block(nov*nbj, nov, datas);
        const int m = nbasis2;
        const int n = jsize * bsize;
        blas::transpose(datas, m, n, data);

        for (int nka = -K_; nka < maxK1; ++nka, ++nbja) {
          const int nkac = nka + K_;
          const std::complex<double> exponent(0.0, K_ != 0 ? (m1 * nka * pi)/ K_ : 0.0);
          const std::complex<double> prefac = exp(exponent);
          const int nsize = jsize * bsize;
//        #pragma omp parallel for
          for (int ii = 0; ii < nbasis1; ++ii) {
            const int offset1 = nsize * nbasis1 * ii;
            const int offset2 = nsize * asize * ii;
            zgemm3m_("N", "N", &nsize, &asize, &nbasis1, &prefac, data + offset1, &nsize,
                                                                coeff_a->bp(nka) + nbasis1 * astart, &nbasis1, &czero,
                                                                datas + offset2, &nsize);
          }
          if (q1_front)
            intermediate_KKK.append(novv, datas);
          else
            intermediate_KKK.add_block(novv * nbja, novv, datas);
        }
      }
    } // end of contraction a for given m1

    if (q1_front)
      intermediate_KKK.reopen_with_inout();

  } // end of m1 loop

  // now intermediate_KKK is ready.
  for (int nkb = -K_; nkb != maxK1; ++nkb) {
    int nbja = (nkb + K_) * KK * KK;
    for (int nkj = -K_; nkj != maxK1; ++nkj) {
      for (int nka = -K_; nka < maxK1; ++nka, ++nbja) {
        intermediate_KKK.get_block(novv * nbja, novv, intermediate_novv);
        std::complex<double>* conjc2 = new std::complex<double>[nbasis1 * isize];

        // momentum conservation
        int nki = nka + nkb - nkj;
        if (nki < - K_) nki += K_ * 2;
        else if (nki >=  K_) nki -= K_ * 2;

        const std::complex<double>* idata = coeff_i->bp(nki) + nbasis1 * istart;
        for (int ii = 0; ii != nbasis1 * isize; ++ii) conjc2[ii] = conj(idata[ii]);
        const int nsize = jsize * asize * bsize;
        zgemm3m_("N", "N", &nsize, &isize, &nbasis1, &cone, intermediate_novv, &nsize,
                                                          conjc2, &nbasis1, &czero,
                                                          datas, &nsize);
        delete[] conjc2;
        mo_int->append(noovv, datas);
      }
    }
  } // end of contraction i

  delete[] data;
  delete[] datas;
  delete[] data_read;
  delete[] conjc;
  delete[] intermediate_novv;
  delete[] blocks;

  std::cout << "  done" << std::endl << std::endl;
  mo_int->reopen_with_inout();
  mo_int->sort_inside_blocks();
  return mo_int;
};

}

#endif

