//
// BAGEL - Parallel electron correlation program.
// Filename: dfinttask_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#ifndef __SRC_DF_DFINTTASK_LONDON_H
#define __SRC_DF_DFINTTASK_LONDON_H

#include <src/df/df_london.h>
#include <src/df/dfblock_london.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/molecule/shell.h>

namespace bagel {

// template <typename TBatch, int N>
// based on a template permitting different integral classes and number of blocks, N
class DFIntTask_London {
  protected:
    constexpr static int N = 1;
    const std::array<std::shared_ptr<const Shell>,4> shell_;
    const std::array<int,3> offset_; // at most 3 elements
    std::array<std::shared_ptr<DFBlock_London>,N> dfblocks_;

    std::shared_ptr<ComplexERIBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<ComplexERIBatch>(input, 2.0);
      eribatch->compute();
      return eribatch;
    }

  public:
    DFIntTask_London(std::array<std::shared_ptr<const Shell>,4>&& a, std::array<int,3>&& b, std::array<std::shared_ptr<DFBlock_London>,N>& df)
     : shell_(a), offset_(b), dfblocks_(df) { };

    void compute();
};


class DFDist_London; // eww

// Reference class could take either DFDist or DFDistT; this one is only needed for DFDist_London, so not templated
// template <typename T>
class DFIntTask_OLD_London {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,2> offset_; // at most 3 elements
    int rank_;
    DFDist_London* df_;

  public:
    DFIntTask_OLD_London(std::array<std::shared_ptr<const Shell>,4>&& a, std::array<int,2>&& b, DFDist_London* df)
     : shell_(a), offset_(b), rank_(offset_.size()), df_(df) { }

    void compute();

};

}

#endif
