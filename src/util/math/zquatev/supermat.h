//
// ZQUATEV: Diagonalization of quaternionic matrices
// File   : supermat.h
// Copyright (c) 2016, Toru Shiozaki (shiozaki@northwestern.edu)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.
//


#ifndef __SUPERMATRIX_H
#define __SUPERMATRIX_H

#include "f77.h"
#include <cassert>
#include <complex>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>

namespace ts {
namespace impl {

template <size_t NB, size_t MB>
class SuperMatrix {
  protected:
    // raw pointer. Efficiency is everything!
    std::complex<double>* data_;

    // size information
    const int nmax_;
    const int mmax_;
    int area() const { return nmax_*mmax_; }

    // size pointer
    std::array<int, NB> nptr_;
    std::array<int, MB> mptr_;

  public:
    SuperMatrix(std::complex<double>* d, const int nm, const int mm, const int nst = 1, const int mst = 1, const bool init = true, const bool last_is_one = false)
        : data_(d), nmax_(nm), mmax_(mm) {
      std::fill(nptr_.begin(), nptr_.end(), nst);
      std::fill(mptr_.begin(), mptr_.end(), mst);
      if (last_is_one)
        mptr_[mst-1] = 1;
      if (init)
        std::fill_n(data_, (!last_is_one ? nmax_*mmax_*NB*MB : (nmax_*mmax_*NB*(MB-1) + nmax_*NB)), 0.0);
    }

    SuperMatrix(std::complex<double>* d, const int nm, const int mm, const std::array<int, NB>& nptr, const std::array<int, MB>& mptr)
        : data_(d), nmax_(nm), mmax_(mm), nptr_(nptr), mptr_(mptr) {
    }

    SuperMatrix(std::complex<double>* d, const SuperMatrix<NB,MB>& o)
        : data_(d), nmax_(o.nmax_), mmax_(o.mmax_), nptr_(o.nptr_), mptr_(o.mptr_) {
      for (int m = 0; m != MB; ++m)
        for (int n = 0; n != NB; ++n)
          for (int j = 0; j != mptr_[m]; ++j)
            std::copy_n(o.block(n,m)+j*nmax_, nptr_[n], block(n,m)+j*nmax_);
    }

    template<int iblock, int mb = 1, class = typename std::enable_if<(iblock+mb-1 < MB)>::type>
    SuperMatrix<NB,mb> slice() {
      std::array<int, mb> mp;
      std::copy_n(mptr_.begin()+iblock, mb, mp.begin());
      return SuperMatrix<NB,mb>(block(0,iblock), nmax_, mmax_, nptr_, mp);
    }

    template<int nb, class = typename std::enable_if<(nb < NB)>::type>
    SuperMatrix<nb,MB> trunc() {
      static_assert(MB == 1, "only defined for MB = 1");
      std::array<int, nb> np;
      std::copy_n(nptr_.begin(), nb, np.begin());
      return SuperMatrix<nb,MB>(block(0,0), nmax_, mmax_, np, mptr_);
    }

    SuperMatrix<NB,MB> shift(const int n) {
      std::array<int, NB> np;
      for (int i = 0; i != NB; ++i) {
        assert(nptr_[i]-n > 0);
        np[i] = nptr_[i]-n;
      }
      return SuperMatrix<NB,MB>(data_+n, nmax_, mmax_, np, mptr_);
    }

    template<size_t I, size_t J, class = typename std::enable_if<(I < NB and J < MB)>::type>
    std::complex<double>& data(const int n, const int m) { return *(block(I, J) + n + nmax_*m); }
    template<size_t I, size_t J, class = typename std::enable_if<(I < NB and J < MB)>::type>
    std::complex<double>* ptr(const int n, const int m) { return block(I, J) + n + nmax_*m; }

    std::complex<double>* block(const int i, const int j) { return data_+area()*(i+NB*j); }
    const std::complex<double>* block(const int i, const int j) const { return data_+area()*(i+NB*j); }

    void conj() {
      for (int m = 0; m != MB; ++m)
        for (int n = 0; n != NB; ++n)
          for (int j = 0; j != mptr_[m]; ++j)
            conj_n(block(n,m)+j*nmax_, nptr_[n]);
    }

    void scale(const std::complex<double>& a) {
      for (int m = 0; m != MB; ++m)
        for (int n = 0; n != NB; ++n)
          for (int j = 0; j != mptr_[m]; ++j)
            zscal_(nptr_[n], a, block(n,m)+j*nmax_, 1);
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < MB)>::type>
    void append_column(const std::complex<double>* d, const int ld, const int off = 0) {
      assert(mptr_[iblock]+1 <= mmax_);
      for (int nblock = 0; nblock != NB; ++nblock)
        std::copy_n(d+ld*nblock, nptr_[nblock]-off, block(nblock, iblock)+mptr_[iblock]*nmax_+off);
      mptr_[iblock]++;
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < MB)>::type>
    void append_column(const SuperMatrix<NB,1>& d) {
      append_column<iblock>(d.block(0,0), d.nmax());
    }

    template<size_t iblock, size_t nblock = 0, class = typename std::enable_if<(iblock < MB)>::type>
    void append_column(const int off = 0, const std::complex<double>& a = 0.0) {
      assert(mptr_[iblock]+1 <= nmax_);
      *(block(nblock, iblock)+off+nmax_*mptr_[iblock]) = a;
      mptr_[iblock]++;
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < MB)>::type>
    void add_lastcolumn(const std::complex<double>* d, const int ld, const int off = 0, const std::complex<double>& a = 1.0) {
      for (int nblock = 0; nblock != NB; ++nblock)
        zaxpy_(nptr_[nblock]-off, a, d+ld*nblock, 1, block(nblock, iblock)+(mptr_[iblock]-1)*nmax_+off, 1);
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < MB)>::type>
    void add_lastcolumn(const SuperMatrix<NB,1>& d, const std::complex<double>& a = 1.0) {
      assert(d.mptr(0) == 1);
      add_lastcolumn<iblock>(d.block(0,0), d.nmax(), 0, a);
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < MB)>::type>
    void write_lastcolumn(std::complex<double>* d, const int ld, const int off = 0) {
      for (int nblock = 0; nblock != NB; ++nblock)
        std::copy_n(d+ld*nblock, nptr_[nblock]-off, block(nblock, iblock)+(mptr_[iblock]-1)*nmax_+off);
    }

    template<size_t iblock, class = typename std::enable_if<(iblock < NB)>::type>
    void append_row(std::complex<double>* d, const int ld) {
      assert(nptr_[iblock]+1 <= nmax_);
      for (int mblock = 0; mblock != NB; ++mblock)
        for (int m = 0; m != mptr_[mblock]; ++m)
          *(block(iblock, mblock)+m*nmax_+nptr_[iblock]) = *(d+m+ld*mblock);
      nptr_[iblock]++;
    }

    template<size_t iblock, size_t mblock = 0, class = typename std::enable_if<(iblock < NB)>::type>
    void append_row(const int off = 0, const std::complex<double>& a = 0.0) {
      assert(nptr_[iblock]+1 <= nmax_);
      *(block(iblock, mblock)+off*nmax_+nptr_[iblock]) = a;
      nptr_[iblock]++;
    }

    template<size_t nblock>
    void cut_row(const int off, SuperMatrix<MB,1>& result) const {
      for (int i = 0; i != MB; ++i) {
        std::complex<double>* out = result.block(i,0);
        result.nptr(i) = mptr(i);
        for (int m = 0; m != mptr_[i]; ++m)
          *out++ = std::conj(*(block(nblock, i)+m*nmax_+off));
      }
    }

    int& nptr(const int i) { return nptr_[i]; }
    int& mptr(const int i) { return mptr_[i]; }
    const int& nptr(const int i) const { return nptr_[i]; }
    const int& mptr(const int i) const { return mptr_[i]; }
    int nmax() const { return nmax_; }
    int mmax() const { return mmax_; }

    void reset() {
      std::fill(nptr_.begin(), nptr_.end(), 1);
      std::fill(mptr_.begin(), mptr_.end(), 1);
      std::fill_n(data_, nmax_*mmax_*NB*MB, 0.0);
    }

    void print() {
      std::cout << std::setprecision(4);
      for (int n = 0; n != NB; ++n)
        for (int m = 0; m != MB; ++m) {
          std::cout << n << " " << m << ":" << std::endl;
          for (int i = 0; i != nptr_[n]; ++i) {
            for (int j = 0; j != mptr_[m]; ++j)
              std::cout << *(block(n, m)+i+nmax_*j);
            std::cout << std::endl;
          }
        }
    }
};


namespace {

enum Trans { _N, _T, _C};

template<Trans TA, Trans TB, size_t N, size_t M, size_t K, size_t L, size_t X, size_t Y>
void contract(std::complex<double> a, const SuperMatrix<N,M>& A, const SuperMatrix<K,L>& B, SuperMatrix<X,Y>& C) {
  const constexpr char* c1 = TA == _N ? "N" : (TA == _T ? "T" : "C");
  const constexpr char* c2 = TB == _N ? "N" : (TB == _T ? "T" : "C");
  const constexpr bool transA = !(TA == _N);
  const constexpr bool transB = !(TB == _N);
  const constexpr int loopblock = transA ? N : M;
  static_assert((transA ? M : N) == X, "A dim wrong");
  static_assert((transB ? K : L) == Y, "B dim wrong");
  static_assert((transA ? N : M) == (transB ? L : K), "AB dim wrong");
  for (int y = 0; y != Y; ++y)
    for (int x = 0; x != X; ++x) {
      for (int l = 0; l != loopblock; ++l) {
        assert((transA ? A.nptr(l) : A.mptr(l)) == (transB ? B.mptr(l) : B.nptr(l)));
        zgemm3m_(c1, c2, (transA ? A.mptr(x) : A.nptr(x)), (transB ? B.nptr(y) : B.mptr(y)), (transA ? A.nptr(l) : A.mptr(l)),
                 a, (transA ? A.block(l, x) : A.block(x, l)), A.nmax(), (transB ? B.block(y, l) : B.block(l, y)), B.nmax(),
                 1.0, C.block(x,y), C.nmax());
      }
      C.nptr(x) = (transA ? A.mptr(x) : A.nptr(x));
      C.mptr(y) = (transB ? B.nptr(y) : B.mptr(y));
      assert(C.nptr(x) <= C.nmax());
      assert(C.mptr(y) <= C.mmax());
    }
}

template<Trans TA, size_t N, size_t M, size_t K, size_t L, size_t X, size_t Y>
void contract(std::complex<double> a, const SuperMatrix<N,M>& A, const SuperMatrix<K,L>& B, SuperMatrix<X,Y>& C) {
  const constexpr char* c1 = TA == _N ? "N" : (TA == _T ? "T" : "C");
  const constexpr bool transA = !(TA == _N);
  const constexpr int loopblock = transA ? N : M;
  static_assert((transA ? M : N) == X, "A dim wrong");
  static_assert(L == Y, "B dim wrong");
  static_assert((transA ? N : M) == K, "AB dim wrong");
  assert(B.mmax() == 1 && C.mmax() == 1);
  for (int y = 0; y != Y; ++y)
    for (int x = 0; x != X; ++x) {
      for (int l = 0; l != loopblock; ++l) {
        assert((transA ? A.nptr(l) : A.mptr(l)) == B.nptr(l));
        zgemv_(c1, A.nptr(transA ? l : x) , A.mptr(transA ? x : l), a, (transA ? A.block(l, x) : A.block(x, l)), A.nmax(), B.block(l,y), 1, 1.0, C.block(x,y), 1);
      }
      C.nptr(x) = (transA ? A.mptr(x) : A.nptr(x));
      C.mptr(y) = B.mptr(y);
      assert(C.nptr(x) <= C.nmax());
      assert(C.mptr(y) <= C.mmax());
    }
}


template<Trans TA, size_t N, size_t M, size_t K, size_t L, size_t X, size_t Y>
void contract_tr(std::complex<double> a, const SuperMatrix<N,M>& A, const SuperMatrix<K,L>& B, SuperMatrix<X,Y>& C, std::complex<double>* work = nullptr) {
#if 1
  const constexpr char* c1 = TA == _N ? "N" : (TA == _T ? "T" : "C");
  const constexpr bool transA = !(TA == _N);
  const constexpr int loopblock = transA ? N : M;
  static_assert((transA ? M : N) == X, "A dim wrong");
  static_assert(L == Y, "B dim wrong");
  static_assert((transA ? N : M) == K, "AB dim wrong");
  assert(B.mmax() == 1 && C.mmax() == 1);
  assert(loopblock == 1 || work);
  for (int y = 0; y != Y; ++y)
    for (int x = 0; x != X; ++x) {
      for (int l = 0; l != loopblock; ++l) {
        assert((transA ? A.nptr(l) : A.mptr(l)) == B.nptr(l));
        if (A.nptr(transA ? l : x) == A.mptr(transA ? x : l)) {
          std::copy_n(B.block(l,y), B.nptr(l), work);
          ztrmv_("U", c1, "N", B.nptr(l), (transA ? A.block(l, x) : A.block(x, l)), A.nmax(), work, 1);
          zaxpy_(B.nptr(l), a, work, 1, C.block(x,y), 1);
        } else {
          zgemv_(c1, A.nptr(transA ? l : x) , A.mptr(transA ? x : l), a, (transA ? A.block(l, x) : A.block(x, l)), A.nmax(), B.block(l,y), 1, 1.0, C.block(x,y), 1);
        }
      }
      C.nptr(x) = (transA ? A.mptr(x) : A.nptr(x));
      C.mptr(y) = B.mptr(y);
      assert(C.nptr(x) <= C.nmax());
      assert(C.mptr(y) <= C.mmax());
    }
#else
   contract<TA>(a, A, B, C);
#endif
}

}
}}

#endif
