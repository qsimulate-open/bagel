//
// BAGEL - Parallel electron correlation program.
// Filename: asd.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __ASD_ASD_H
#define __ASD_ASD_H

#include <src/asd/asd_base.h>
#include <src/dimer/dimer_prop.h>

namespace bagel {

/// Template for ASD
template <class VecType>
class ASD : public ASD_base {
  protected: using DSubSpace = DimerSubspace<VecType>;
  protected: using DCISpace = DimerCISpace_base<VecType>;
  protected: using CiType = typename VecType::Ci;

  protected:
    std::shared_ptr<DCISpace> cispace_;
    std::vector<DSubSpace> subspaces_;

    std::vector<DimerSubspace_base> subspaces_base() const override {
      std::vector<DimerSubspace_base> out; out.reserve(subspaces_.size());
      for (auto& i : subspaces_) out.push_back(DimerSubspace_base(i));
      return out;
    }

  public:
    ASD(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DCISpace> cispace);

    void compute() override;

  private:

    std::shared_ptr<Matrix> compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const;
    std::shared_ptr<Matrix> compute_diagonal_1e(const DSubSpace& subspace, const double* hAA, const double* hBB, const double diag) const;

    // Diagonal block stuff
    void compute_pure_terms(DSubSpace& subspace, std::shared_ptr<const DimerJop> jop);
    void compute_diagonal_2rdm(DSubSpace& subspace);

    virtual std::shared_ptr<VecType> form_sigma(std::shared_ptr<const VecType> ccvec, std::shared_ptr<const MOFile> jop) const = 0;
    virtual std::shared_ptr<VecType> form_sigma_1e(std::shared_ptr<const VecType> ccvec, const double* modata) const = 0;

};

// Locks to make sure the following files are not included on their own
#define ASD_HEADERS
#include <src/asd/asd_compute.hpp>
#include <src/asd/asd_compute_diagonal.hpp>
#include <src/asd/asd_init.hpp>
#undef ASD_HEADERS

}

#endif
