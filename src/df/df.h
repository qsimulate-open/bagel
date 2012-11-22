//
// BAGEL - Parallel electron correlation program.
// Filename: df.h
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

#ifndef __SRC_DF_DF_H
#define __SRC_DF_DF_H

#include <vector>
#include <memory>
#include <list>
#include <map>
#include <stdexcept>
#include <stddef.h>
#include <src/df/dfblock.h>
#include <src/scf/atom.h>
#include <src/util/matrix.h>
#include <src/df/dfinttask_old.h>

namespace bagel {

class DFHalfDist;
class DFFullDist;

class ParallelDF {
  protected:
    // blocks that this process has
    std::shared_ptr<DFBlock> block_;
    // hash key and process number
    std::map<int, int> global_table_;
    std::vector<std::pair<int, int> > atable_;
    size_t maxasize_;

  public:
    ParallelDF();

    void add_block(std::shared_ptr<DFBlock> o);

    std::unique_ptr<double[]> form_2index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::unique_ptr<double[]> form_4index(std::shared_ptr<const ParallelDF> o, const double a, const bool swap = false) const;
    std::shared_ptr<Matrix> form_aux_2index(std::shared_ptr<const ParallelDF> o, const double a) const; 

    void daxpy(const double a, const std::shared_ptr<const ParallelDF> o);
    void scale(const double a);

    std::unique_ptr<double[]> get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const;

    size_t maxasize() const { return maxasize_; }
    const std::pair<int, int> atable(const int i) const { return atable_[i]; }

};


class DFDist : public ParallelDF, public std::enable_shared_from_this<DFDist> {
  friend class DFIntTask_OLD<DFDist>;
  friend class DFFullDist;
  friend class DFHalfDist;
  protected:
    // #orbital basis
    const size_t nbasis0_; // outer
    const size_t nbasis1_; // inner
    // #auxiliary basis
    const size_t naux_;

    // AO two-index integrals ^ -1/2
    std::shared_ptr<Matrix> data2_;

    std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input);

    void common_init(const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&,
                     const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv);
    void make_table(const int nmax);

  public:
    // construction of a block from AO integrals
    DFDist(const int nbas, const int naux, const std::vector<std::shared_ptr<const Atom> >& atoms,
                                           const std::vector<std::shared_ptr<const Atom> >& aux_atoms, const double thr, const bool inverse, const double dum)
      : nbasis0_(nbas), nbasis1_(nbas), naux_(naux)  {
      common_init(atoms, atoms, aux_atoms, thr, inverse);
    }

    DFDist(const std::shared_ptr<const DFDist> df) : nbasis0_(df->nbasis0_), nbasis1_(df->nbasis1_), naux_(df->naux_) { }

    // This might not be beautiful - for historical reasons. Will remove later TODO
    DFDist(const int nbas0, const int nbas1, const int naux, const std::vector<const double*> cd, const std::vector<const double*> dd);

    bool has_2index() const { return data2_.get() != nullptr; };
    size_t nbasis0() const { return nbasis0_; };
    size_t nbasis1() const { return nbasis1_; };
    size_t naux() const { return naux_; };

    void add_direct_product(std::vector<const double*> a, std::vector<const double*> b, const double fac);

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DFHalfDist> compute_half_transform(const double* c, const size_t nocc) const;

    // compute a J operator, given density matrices in AO basis
    std::unique_ptr<double[]> compute_Jop(const double* den) const;

    std::unique_ptr<double[]> compute_cd(const double* den) const;

};


class DFHalfDist : public ParallelDF {
  friend class DFFullDist;
  protected:
    const std::shared_ptr<const DFDist> df_;
    const size_t nocc_; // inner
    const size_t nbasis_; // outer
    const size_t naux_;

    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const Matrix> o) const;

  public:
    DFHalfDist(const std::shared_ptr<const DFDist> df, const int nocc) : df_(df), nocc_(nocc), nbasis_(df_->nbasis0()), naux_(df_->naux()) { }

    size_t nocc() const { return nocc_; };
    size_t nbasis() const { return nbasis_; };
    size_t naux() const { return naux_; };
    size_t size() const { return naux_*nbasis_*nocc_; };

    std::shared_ptr<DFFullDist> compute_second_transform(const double* c, const size_t nocc) const;
    std::shared_ptr<DFDist> back_transform(const double* c) const;

    std::shared_ptr<DFHalfDist> copy() const; 
    std::shared_ptr<DFHalfDist> clone() const; 

    void rotate_occ(const double* d);
    std::shared_ptr<DFHalfDist> apply_density(const double* d) const;

    std::unique_ptr<double[]> compute_Kop_1occ(const double* den) const;

    std::shared_ptr<DFHalfDist> apply_J() const { return apply_J(df_->data2_); }
    std::shared_ptr<DFHalfDist> apply_JJ() const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*df_->data2_**df_->data2_))); }
    std::shared_ptr<DFHalfDist> apply_J(const std::shared_ptr<const DFDist> d) const { return apply_J(d->data2_); }
    std::shared_ptr<DFHalfDist> apply_JJ(const std::shared_ptr<const DFDist> d) const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*d->data2_**d->data2_))); }

    const std::shared_ptr<const DFDist> df() const { return df_; }
};


class DFFullDist : public ParallelDF {
  protected:
    const std::shared_ptr<const DFDist> df_;
    const size_t nocc1_; // inner
    const size_t nocc2_; // outer
    const size_t naux_;

    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const Matrix> o) const;

  public:
    DFFullDist(const std::shared_ptr<const DFDist> df, const int nocc1, const int nocc2) : df_(df), nocc1_(nocc1), nocc2_(nocc2), naux_(df_->naux()) { }

    int nocc1() const { return nocc1_; }
    int nocc2() const { return nocc2_; }
    int naux() const { return naux_; }

    std::shared_ptr<DFFullDist> copy() const; 
    std::shared_ptr<DFFullDist> clone() const; 

    std::shared_ptr<DFHalfDist> back_transform(const double* c) const;

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

    std::unique_ptr<double[]> form_4index_1fixed(const std::shared_ptr<const DFFullDist> o, const double a, const size_t n) const;

    // utility functions
    std::shared_ptr<Matrix> form_aux_2index_apply_J(const std::shared_ptr<const DFFullDist> o, const double a) const;
    void set_product(const std::shared_ptr<const DFFullDist>, const std::unique_ptr<double[]>&, const int jdim, const size_t offset); 

    std::shared_ptr<DFFullDist> apply_J() const { return apply_J(df_->data2_); }
    std::shared_ptr<DFFullDist> apply_JJ() const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*df_->data2_**df_->data2_))); }
    std::shared_ptr<DFFullDist> apply_J(const std::shared_ptr<const DFDist> d) const { return apply_J(d->data2_); }
    std::shared_ptr<DFFullDist> apply_JJ(const std::shared_ptr<const DFDist> d) const { return apply_J(std::shared_ptr<Matrix>(new Matrix(*d->data2_**d->data2_))); }

    const std::shared_ptr<const DFDist> df() const { return df_; }
};

}


#endif
