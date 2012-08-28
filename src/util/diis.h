//
// BAGEL - Parallel electron correlation program.
// Filename: diis.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __src_util_diis_h
#define __src_util_diis_h

#include <list>
#include <src/util/f77.h>
#include <cassert>
#include <memory>

// std::shared_ptr<T> is assumed to be a shared_pointer of some class
// which have daxpy and ddot functions.
// T must have clone() function that returns shared_ptr<T>

template <class T>
class DIIS {
  typedef std::shared_ptr<const T> RefT;
  typedef std::list<std::pair<RefT, RefT> > Container_type_;
  typedef typename Container_type_::iterator iterator;

  protected:
    const int ndiis_;
    const int nld_;
    const int lwork_;

    Container_type_ data_;

    std::unique_ptr<double[]> matrix_;
    std::unique_ptr<double[]> matrix_save_;
    std::unique_ptr<double[]> coeff_;
    std::unique_ptr<double[]> work_;
    std::unique_ptr<int[]> ipiv_;

  public:
    DIIS(const int ndiis) : ndiis_(ndiis), nld_(ndiis+1), lwork_((ndiis+1)*(ndiis+1)),
      matrix_(new double[(ndiis+1)*(ndiis+1)]),
      matrix_save_(new double[(ndiis+1)*(ndiis+1)]),
      coeff_(new double[ndiis+1]),
      work_(new double[(ndiis+1)*(ndiis+1)]),
      ipiv_(new int[(ndiis+1)*(ndiis+1)]) { };

    ~DIIS() { };

    std::shared_ptr<T> extrapolate(const std::pair<RefT, RefT> input) {
      RefT v = input.first;
      RefT e = input.second;
      data_.push_back(input);
      const int size = nld_ * nld_;

      if (data_.size() > ndiis_) {
        data_.pop_front();
        for (int i = 1; i != ndiis_; ++i) {
          for (int j = 1; j != ndiis_; ++j)
            matrix_[(j - 1) + (i - 1) * nld_] = matrix_save_[j + i * nld_];
        }
      } else if (data_.size() != 1) {
        dcopy_(size, matrix_save_, 1, matrix_, 1);
      }
      const int cnum = data_.size();
      iterator data_iter = data_.begin();

      for (int i = 0; i != cnum - 1; ++i, ++data_iter)
        matrix_[(cnum - 1) + i * nld_] = matrix_[i + (cnum - 1) * nld_] = e->ddot(*(data_iter->second));
      matrix_[(cnum - 1) + (cnum - 1) * nld_] = e->ddot(e);
      for (int i = 0; i != cnum; ++i)
        matrix_[cnum + i * nld_] = matrix_[i + cnum * nld_] = -1.0;
      matrix_[cnum + cnum * nld_] = 0.0;
      for (int i = 0; i != cnum; ++i) coeff_[i] = 0.0;
      coeff_[cnum] = -1.0;

      dcopy_(size, matrix_, 1, matrix_save_, 1);

      const int cdim = cnum + 1;
      int info;
      dsysv_("U", cdim, 1, matrix_, nld_, ipiv_, coeff_, nld_, work_, lwork_, info);
      if (info) throw std::runtime_error("DSYSV failed in diis.h");

      std::shared_ptr<T> out = input.first->clone();

      data_iter = data_.begin();
      for (int i = 0; i != cnum; ++i, ++data_iter)
        out->daxpy(coeff_[i], *(data_iter->first));

      return out;
    };

};

#endif

