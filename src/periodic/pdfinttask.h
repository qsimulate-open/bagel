//
// BAGEL - Parallel electron correlation program.
// Filename: pdfinttask.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#ifndef __SRC_PERIODIC_PDFINTTASK_H
#define __SRC_PERIODIC_PDFINTTASK_H

#include <src/df/dfblock.h>
#include <src/molecule/shell.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/naibatch.h>
#include <src/integral/os/overlapbatch.h>

namespace bagel {

class PDFIntTask_aux {
  protected:
    std::array<std::shared_ptr<const Shell>,2> shell_;
    int offset_;
    std::shared_ptr<VectorB> pdf_;

  public:
    // <i|.>
    PDFIntTask_aux(std::array<std::shared_ptr<const Shell>,2>&& sh, int offset, std::shared_ptr<VectorB> pdf)
     : shell_(sh), offset_(offset), pdf_(pdf) { }

    void compute() {

      auto overlap = std::make_shared<OverlapBatch>(shell_);
      const double* odata = overlap->data();
      double* const data = pdf_->data();

      for (int j0 = offset_; j0 != offset_ + shell_[1]->nbasis(); ++j0, ++odata)
        data[j0] = *odata;
    }
};


class PDFIntTask_coeff {
  protected:
    std::array<std::shared_ptr<const Shell>,2> shell_;
    std::array<int,2> offset_;
    std::shared_ptr<btas::Tensor3<double>> coeffC_;
    std::shared_ptr<const VectorB> data1_;
    int ncell_;

    std::shared_ptr<OverlapBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,2>& input) const {
      auto obatch = std::make_shared<OverlapBatch>(input);
      obatch->compute();
      return obatch;
    }

  public:
    // <r|sL>
    PDFIntTask_coeff(std::array<std::shared_ptr<const Shell>,2>&& sh, std::array<int,2>&& offset,
                     std::shared_ptr<btas::Tensor3<double>>& coeffC, const std::shared_ptr<const VectorB>& data1, const int n)
     : shell_(sh), offset_(offset), coeffC_(coeffC), data1_(data1), ncell_(n) { }

    void compute() {

      std::shared_ptr<OverlapBatch> overlap = compute_batch(shell_);
      const double* odata = overlap->data();

      const size_t naux = coeffC_->extent(0);
      const size_t nbas = coeffC_->extent(1);
      const double q = data1_->rms() * naux;

      double* const coeff = coeffC_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[1]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++odata)
          for (int a = 0; a != naux; ++a)
            coeff[a + naux * (j1 + nbas * j0)] = *odata * (*data1_)(a) * ncell_ / q;

    }
};


class PDFIntTask_2index {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,2> offset_;
    std::shared_ptr<Matrix> data2_;

    std::shared_ptr<ERIBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<ERIBatch>(input, 2.0);
      eribatch->compute();
      return eribatch;
    }

  public:
    // (i.|jL.)
    PDFIntTask_2index(std::array<std::shared_ptr<const Shell>,4>&& sh, std::array<int,2>&& offset, std::shared_ptr<Matrix>& data2)
     : shell_(sh), offset_(offset), data2_(data2) { }

    void compute() {

      std::shared_ptr<ERIBatch> eribatch = compute_batch(shell_);
      const double* eridata = eribatch->data();

      const size_t naux = data2_->ndim();

      assert(offset_.size() == 2);
      double* const data = data2_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++eridata)
          data[j1+j0*naux] = *eridata;
    }
};


class PDFIntTask_3index {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,3> offset_;
    std::shared_ptr<DFBlock> dfblock_;

    std::shared_ptr<ERIBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<ERIBatch>(input, 2.0);
      eribatch->compute();
      return eribatch;
    }

  public:
    // (r sL'|iL .)
    PDFIntTask_3index(std::array<std::shared_ptr<const Shell>,4>&& shells, std::array<int,3>&& offset, std::shared_ptr<DFBlock>& df)
     : shell_(shells), offset_(offset), dfblock_(df) { };

    void compute() {

      std::shared_ptr<ERIBatch> eribatch = compute_batch(shell_);

      assert(dfblock_->b1size() == dfblock_->b2size());
      assert(offset_.size() == 3);
      const size_t nbin = dfblock_->b1size();
      const size_t naux = dfblock_->asize();
      const double* eridata = eribatch->data();

      double* const data = dfblock_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1, eridata += shell_[1]->nbasis())
          std::copy_n(eridata, shell_[1]->nbasis(), data + offset_[2] + naux * (j1 + nbin * j0));
    }
};


class PDFIntTask_NAI {

  protected:
    // (r sL'|\delta_L)
    std::array<std::shared_ptr<const Shell>,2> shell_;
    std::array<int,2> offset_;
    std::shared_ptr<const Molecule> mol_;
    std::shared_ptr<Matrix> data2_;

    std::shared_ptr<NAIBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,2>& input,
                                            const std::shared_ptr<const Molecule> mol) const {
      auto naibatch = std::make_shared<NAIBatch>(input, mol);
      naibatch->compute();
      return naibatch;
    }

  public:
    PDFIntTask_NAI(std::array<std::shared_ptr<const Shell>,2>&& shells, std::array<int,2>&& offset,
                   const std::shared_ptr<const Geometry> mol, std::shared_ptr<Matrix>& data2)
     : shell_(shells), offset_(offset), mol_(mol), data2_(data2) { };

    void compute() {

      std::shared_ptr<NAIBatch> naibatch = compute_batch(shell_, mol_);
      const double* naidata = naibatch->data();
      const int nbas = data2_->ndim();

      assert(offset_.size() == 2);
      double* const data = data2_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[1]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++naidata)
          data[j1+j0*nbas] = *naidata;
    }
};

}

#endif
