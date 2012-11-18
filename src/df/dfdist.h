//
// BAGEL - Parallel electron correlation program.
// Filename: dfdist.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#ifndef __SRC_DF_DFDIST_H
#define __SRC_DF_DFDIST_H

#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <src/rysint/libint.h>

namespace bagel {

class DFHalfDist;
class DFFullDist;


class ParallelDF {
  protected:
    std::list<std::shared_ptr<DFBlock> > blocks_;
  public:
    ParallelDF() {};

    void add_block(std::shared_ptr<DFBlock> o) { blocks_.push_back(o); }

    std::unique_ptr<double[]> form_2index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::unique_ptr<double[]> form_4index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;

};


class DFDist : public DensityFit, public ParallelDF {
  protected:
    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override;

    void common_init(const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv) override;

  public:
    // construction of a block from AO integrals
    DFDist(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom> >& atoms,
                                           const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thr, const bool inverse)
      : DensityFit(nbas, naux) {
      common_init(atoms, atoms, aux_atoms, thr, inverse);
      assert(data_ == nullptr);
    }

    DFDist(const std::shared_ptr<const DensityFit> df) : DensityFit(df->nbasis1(), df->naux()) {
      assert(df->nbasis0() == df->nbasis1());
    }

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist> compute_half_transform(const double* c, const size_t nocc) const;

    // compute a J operator, given density matrices in AO basis
    std::unique_ptr<double[]> compute_Jop(const double* den) const override;

    std::unique_ptr<double[]> compute_cd(const double* den) const override;

};


class DFHalfDist : public DF_Half, public ParallelDF {
  friend class DFFullDist;
  protected:

  public:
    DFHalfDist(const std::shared_ptr<const DensityFit> df, const int nocc) : DF_Half(df, nocc) {
      assert(data_ == nullptr);
    }

    std::shared_ptr<DFFullDist> compute_second_transform(const double* c, const size_t nocc) const;
    std::shared_ptr<DFDist> back_transform(const double* c) const;

    std::shared_ptr<DFHalfDist> copy() const; 
    std::shared_ptr<DFHalfDist> clone() const; 

    void rotate_occ(const double* d);
    std::shared_ptr<DFHalfDist> apply_density(const double* d) const;
};


class DFFullDist : public DF_Full, public ParallelDF {
  protected:

  public:
    DFFullDist(const std::shared_ptr<const DensityFit> df, const int nocc0, const int nocc1) : DF_Full(df, nocc0, nocc1) {
      assert(data_ == nullptr);
    }

    std::shared_ptr<DFFullDist> copy() const; 
    std::shared_ptr<DFFullDist> clone() const; 

    std::shared_ptr<DFHalfDist> back_transform(const double* c) const;

    void daxpy(const double a, const DFHalfDist& o);
    void daxpy(const double a, const DFFullDist& o);
    void scale(const double a);

    void symmetrize();

    // 2RDM contractions
    // special function for RHF
    std::shared_ptr<DFFullDist> apply_closed_2RDM() const;
    // special function for UHF
    std::shared_ptr<DFFullDist> apply_uhf_2RDM(const double* rdma, const double* rdmb) const;
    // general case with closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const;
    // general case without closed orbitals
    std::shared_ptr<DFFullDist> apply_2rdm(const double* rdm) const;
};

}


#endif
