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


#ifndef __NEWINT_DF_DensityFit_H
#define __NEWINT_DF_DensityFit_H

#include <vector>
#include <memory>
#include <stddef.h>
#include <src/df/dfblock.h>
#include <src/scf/atom.h>
#include <src/util/f77.h>
#include <src/rysint/eribatch.h>
#include <src/util/matrix.h>
#include <stdexcept>

namespace bagel {

class DF_Half;
class DF_Full;


class DensityFit : public std::enable_shared_from_this<DensityFit> {
  friend class DFIntTask_OLD;
  friend class DF_Half;
  friend class DF_Full;

  protected:
    // #orbital basis
    const size_t nbasis0_; // outer
    const size_t nbasis1_; // inner
    // #auxiliary basis
    const size_t naux_;

    // AO three-index integrals (not set in DFDist).
    std::shared_ptr<DFBlock> data_;
    // AO two-index integrals ^ -1/2
    std::shared_ptr<Matrix> data2_;

    virtual void common_init(const std::vector<std::shared_ptr<const Atom> >&,
                             const std::vector<std::shared_ptr<const Atom> >&,
                             const std::vector<std::shared_ptr<const Atom> >&, const double thresh, const bool compute_inv);

    // returns a pointer to a stack memory area
    virtual std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) = 0;

    size_t size() const { return nbasis0_*nbasis1_*naux_; };

  public:
    DensityFit(const int nbas, const int naux) : nbasis0_(nbas), nbasis1_(nbas), naux_(naux) {};
    DensityFit(const int nbas0, const int nbas1, const int naux) : nbasis0_(nbas0), nbasis1_(nbas1), naux_(naux) {};
    ~DensityFit() {};

    bool has_2index() const { return data2_.get() != nullptr; };

    size_t nbasis0() const { return nbasis0_; };
    size_t nbasis1() const { return nbasis1_; };
    size_t naux() const { return naux_; };

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DF_Half> compute_half_transform(const double* c, const size_t nocc) const;

    // compute a J operator, given density matrices in AO basis
    virtual std::unique_ptr<double[]> compute_Jop(const double* den) const;

    virtual std::unique_ptr<double[]> compute_cd(const double* den) const;

};


class DF_AO : public DensityFit {
  protected:
    virtual std::pair<const double*, std::shared_ptr<RysInt> > compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override {
      throw std::logic_error("DF_AO::compute_batch should not be called");
      return std::pair<const double*, std::shared_ptr<RysInt> >();
    };
  public:
    DF_AO(const int nbas0, const int nbas1, const int naux, std::unique_ptr<double[]>& dat) : DensityFit(nbas0, nbas1, naux) {
      data_ = std::shared_ptr<DFBlock>(new DFBlock(dat, naux, nbas1, nbas0, 0,0,0));
    }
    DF_AO(const int nbas0, const int nbas1, const int naux, std::shared_ptr<DFBlock> dat) : DensityFit(nbas0, nbas1, naux) {
      data_ = dat;
    }
    // contructor for a seperable part of nuclear gradients
    DF_AO(const int nbas0, const int nbas1, const int naux, const std::vector<const double*> cd, const std::vector<const double*> dd);
    ~DF_AO() {};

// TODO this will be removed.
#if 1
    double* ptr(const size_t i, const size_t j, const size_t k) { return data_->get()+i+naux_*(j+nbasis1_*k); };
    const double* ptr(const size_t i, const size_t j, const size_t k) const { return data_->get()+i+naux_*(j+nbasis1_*k); };
#endif

    void daxpy(const double a, const std::shared_ptr<const DF_AO> o) { daxpy_(size(), a, o->data_->get(), 1, data_->get(), 1); };
};


class DF_Half {
  friend class DF_Full;

  protected:
    const std::shared_ptr<const DensityFit> df_;
    const size_t nocc_;

    std::shared_ptr<DFBlock> data_;
    const size_t nbasis_;
    const size_t naux_;

  public:
    DF_Half(const std::shared_ptr<const DensityFit> df, const int nocc, std::unique_ptr<double[]>& in)
     : df_(df), nocc_(nocc), data_(new DFBlock(in, naux_, nocc, nbasis_, 0,0,0)), nbasis_(df->nbasis0()), naux_(df->naux()) {};
    DF_Half(const std::shared_ptr<const DensityFit> df, const int nocc, std::shared_ptr<DFBlock> in)
     : df_(df), nocc_(nocc), data_(in), nbasis_(df->nbasis0()), naux_(df->naux()) {};
    DF_Half(const std::shared_ptr<const DensityFit> df, const int nocc)
     : df_(df), nocc_(nocc), nbasis_(df->nbasis0()), naux_(df->naux()) {};

    ~DF_Half() {};

    size_t nocc() const { return nocc_; };
    size_t nbasis() const { return nbasis_; };
    size_t naux() const { return naux_; };
    size_t size() const { return naux_*nbasis_*nocc_; };
    const std::shared_ptr<const DensityFit> df() { return df_; };

    std::shared_ptr<DF_Half> clone() const;
    std::shared_ptr<DF_Half> copy() const;

    std::shared_ptr<DF_Half> apply_J() const { return apply_J(df_); };
    std::shared_ptr<DF_Half> apply_JJ() const { return apply_JJ(df_); };
    std::shared_ptr<DF_Half> apply_J(std::shared_ptr<const DensityFit> d) const;
    std::shared_ptr<DF_Half> apply_JJ(std::shared_ptr<const DensityFit> d) const;
    std::shared_ptr<DF_Half> apply_density(const double* d) const;

    std::shared_ptr<DF_Full> compute_second_transform(const double* c, const size_t nocc) const;

    std::shared_ptr<DFBlock> data() const { return data_; };

    // form 2 index quantities
    std::unique_ptr<double[]> form_2index(const std::shared_ptr<const DF_Full> o, const double a) const;
    std::unique_ptr<double[]> form_2index(const std::shared_ptr<const DF_Half> o, const double a) const;
    // form 2 index quantities by contracting Naux and Nao (targeting an nocc*nbasis matrix)
    std::unique_ptr<double[]> form_2index(const std::shared_ptr<const DensityFit> o, const double a) const;

    std::shared_ptr<Matrix> form_aux_2index(const std::shared_ptr<const DF_Half> o) const;

    // form K^ij_rs operator
    std::unique_ptr<double[]> form_4index() const;

    // compute a K operator with one occupied index K_rj(D_tu), given an AO density matrix.
    virtual std::unique_ptr<double[]> compute_Kop_1occ(const double* den) const;

    // AO back transformation
    std::shared_ptr<DF_AO> back_transform(const double*) const;

    void rotate_occ(const double*); 

};


class DF_Full {
  friend class DF_Half;

  protected:
    const std::shared_ptr<const DensityFit> df_;
    const size_t nocc1_; // inner
    const size_t nocc2_; // outer

    std::shared_ptr<DFBlock> data_;
    const size_t naux_;

  public:
    DF_Full(const std::shared_ptr<const DensityFit> df, const size_t nocc1, const size_t nocc2, std::unique_ptr<double[]>& in)
      : df_(df), nocc1_(nocc1), nocc2_(nocc2), data_(new DFBlock(in, naux_, nocc1_, nocc2_, 0,0,0)), naux_(df->naux()) {};

    DF_Full(const std::shared_ptr<const DensityFit> df, const size_t nocc1, const size_t nocc2, std::shared_ptr<DFBlock> in)
      : df_(df), nocc1_(nocc1), nocc2_(nocc2), data_(in), naux_(df->naux()) {};

    DF_Full(const std::shared_ptr<const DensityFit> df, const size_t nocc1, const size_t nocc2)
      : df_(df), nocc1_(nocc1), nocc2_(nocc2), naux_(df->naux()) {};

    DF_Full(std::shared_ptr<DF_Half> o) : df_(o->df()), nocc1_(o->nocc()), nocc2_(o->nbasis()), data_(o->data_), naux_(o->naux()) {};

    ~DF_Full() {};

    int nocc1() const { return nocc1_; };
    int nocc2() const { return nocc2_; };
    int naux() const { return naux_; };
    size_t size() const { return nocc1_*nocc2_*naux_; };

    std::shared_ptr<DF_Full> apply_J() const { return apply_J(df_); };
    std::shared_ptr<DF_Full> apply_JJ() const { return apply_JJ(df_); };
    std::shared_ptr<DF_Full> apply_J(std::shared_ptr<const DensityFit> d) const;
    std::shared_ptr<DF_Full> apply_JJ(std::shared_ptr<const DensityFit> d) const;

    // compute Gamma(kl,ij) * (ij|D)
    std::shared_ptr<DF_Full> apply_2rdm(const double* rdm) const;
    std::shared_ptr<DF_Full> apply_2rdm(const double* rdm2, const double* rdm1, const int nclosed, const int nact) const;
    // special function for RHF
    std::shared_ptr<DF_Full> apply_closed_2RDM() const;
    // special function for UHF
    std::shared_ptr<DF_Full> apply_uhf_2RDM(const double* rdma, const double* rdmb) const;

    // forming all internal 4-index MO integrals
    std::unique_ptr<double[]> form_4index() const;
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DF_Full> o) const;

    // the same as above, but with one index fixed at n (of nocc2).
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DF_Full> o, const size_t n) const;

    // forming two-external K operators
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DensityFit> o) const;

    // contract ii and form (D|E)
    std::shared_ptr<Matrix> form_aux_2index(const std::shared_ptr<const DF_Full> o) const;
    std::shared_ptr<Matrix> form_aux_2index_apply_J(const std::shared_ptr<const DF_Full> o) const;

    std::unique_ptr<double[]> form_2index(const std::shared_ptr<const DF_Half> o, const double a);

    std::unique_ptr<double[]> form_2index(const std::shared_ptr<const DF_Full> o, const double a);

    const std::shared_ptr<const DensityFit> df() const { return df_; };

    // compute (D|ia)(ia|j) and set to the location specified by the offset
    void set_product(const std::shared_ptr<const DF_Full>, const std::unique_ptr<double[]>&, const int jdim, const size_t offset); 

    std::shared_ptr<DF_Full> clone() const;
    std::shared_ptr<DF_Full> copy() const;
    void daxpy(const double a, std::shared_ptr<const DF_Full> o) { daxpy(a, *o); };
    void daxpy(const double a, const DF_Full& o);
    void daxpy(const double a, std::shared_ptr<const DF_Half> o) { daxpy(a, *o); };
    void daxpy(const double a, const DF_Half& o);
    DF_Full& operator+=(const DF_Full& o) { daxpy(1.0, o); return *this; };
    DF_Full& operator-=(const DF_Full& o) { daxpy(-1.0, o); return *this; };

    void scale(const double a);
    DF_Full& operator*=(const double a) { scale(a); return *this; };
    DF_Full& operator/=(const double a) { scale(1.0/a); return *this; };

    void symmetrize();

    // AO back transformation for gradient evaluations
    std::shared_ptr<DF_Half> back_transform(const double* c) const;

};

}

#endif

