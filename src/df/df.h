//
// Newint - Parallel electron correlation program.
// Filename: df.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __NEWINT_DF_DensityFit_H
#define __NEWINT_DF_DensityFit_H

#include <vector>
#include <memory>
#include <src/scf/atom.h>
#include <src/util/f77.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>

class DF_Half;
class DF_Full;


class DensityFit : public std::enable_shared_from_this<DensityFit> {
  protected:
    // AO three-index integrals (naux, nbasis1, nbasis0);
    std::unique_ptr<double[]> data_;
    // AO two-index integrals
    std::unique_ptr<double[]> data2_;
    // #orbital basis
    const size_t nbasis0_;
    const size_t nbasis1_;
    // #auxiliary basis
    const size_t naux_;

    // returns three-index integrals
    double* data() { return data_.get(); };
    // returns two-index integrals
    double* data2() { return data2_.get(); };
    const double* const data() const { return data_.get(); };
    const double* const data2() const { return data2_.get(); };

    void common_init(const std::vector<std::shared_ptr<Atom> >&, const std::vector<std::vector<int> >&,
                     const std::vector<std::shared_ptr<Atom> >&, const std::vector<std::vector<int> >&,
                     const std::vector<std::shared_ptr<Atom> >&, const std::vector<std::vector<int> >&, const double, const bool);

    virtual void compute_batch(std::unique_ptr<double[]>&, std::vector<std::shared_ptr<Shell> >&,
                               const int, const int, const int, const int, const int) = 0;

  public:
    DensityFit(const int nbas, const int naux) : nbasis0_(nbas), nbasis1_(nbas), naux_(naux) {};
    ~DensityFit() {};

    const double* const data_3index() const { return data(); };
    const double* const data_2index() const { return data2(); };

    size_t nbasis0() const { return nbasis0_; };
    size_t nbasis1() const { return nbasis1_; };
#if 0
    size_t nbasis() const {
      if (nbasis0_ == nbasis1_) {
        return nbasis0_;
      } else {
        throw std::logic_error("I used to use this function to get #basis, but it should refer to a Geometry object.");
      }
    };
#endif
    size_t naux() const { return naux_; };

    // compute half transforms; c is dimensioned by nbasis_;
    std::shared_ptr<DF_Half> compute_half_transform(const double* c, const size_t nocc) const;

}; 


class DF_Half {

  protected:
    std::unique_ptr<double[]> data_;
    const std::shared_ptr<const DensityFit> df_;
    const size_t nocc_;
    const size_t nbasis_;
    const size_t naux_;

  public:
    DF_Half(const std::shared_ptr<const DensityFit> df, const int nocc, std::unique_ptr<double[]>& in)
     : df_(df), nocc_(nocc), data_(std::move(in)), nbasis_(df->nbasis0()), naux_(df->naux()) {}; 

    ~DF_Half() {};

    double* const data() { return data_.get(); };
    const double* const data() const { return data_.get(); };
    std::unique_ptr<double[]> move_data() { return std::move(data_); };

    std::shared_ptr<DF_Half> clone() const {
      std::unique_ptr<double[]> dat(new double[nocc_*naux_*nbasis_]);
      std::shared_ptr<DF_Half> out(new DF_Half(df_, nocc_, dat)); 
      return out;
    };

    std::shared_ptr<DF_Half> apply_J() const {
      std::unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
      dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, df_->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_); 
      std::shared_ptr<DF_Half> tmp(new DF_Half(df_, nocc_, out));
      return tmp; 
    } 

    std::shared_ptr<DF_Full> compute_second_transform(const double* c, const size_t nocc) const;

    void form_2index(std::unique_ptr<double[]>& target, const double a = 1.0, const double b = 0.0) const;
    void form_2index(std::unique_ptr<double[]>& target, std::shared_ptr<const DF_Full> o, const double a = 1.0, const double b = 0.0) const;

    // form K operator
    std::unique_ptr<double[]> form_4index() const;
    void form_4index(std::unique_ptr<double[]>& target) const;
};


class DF_Full {

  protected:
    std::unique_ptr<double[]> data_;
    const std::shared_ptr<const DensityFit> df_;
    const size_t nocc1_; // inner
    const size_t nocc2_; // outer
    const size_t naux_;

  public:
    DF_Full(const std::shared_ptr<const DensityFit> df, const size_t nocc1, const size_t nocc2, std::unique_ptr<double[]>& in)
      : df_(df), nocc1_(nocc1), nocc2_(nocc2), data_(std::move(in)), naux_(df->naux()) {};

    ~DF_Full() {};

    const int nocc1() const { return nocc1_; };
    const int nocc2() const { return nocc2_; };

    std::shared_ptr<DF_Full> apply_J() const {
      std::unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
      dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, df_->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_); 
      std::shared_ptr<DF_Full> tmp(new DF_Full(df_, nocc1_, nocc2_, out));
      return tmp; 
    };

    std::shared_ptr<DF_Full> apply_2rdm(const double* rdm) {
      std::unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
      dgemm_("N", "T", naux_, nocc1_*nocc2_, nocc1_*nocc2_, 1.0, data_.get(), naux_, rdm, nocc1_*nocc2_, 0.0, out.get(), naux_); 
      std::shared_ptr<DF_Full> tmp(new DF_Full(df_, nocc1_, nocc2_, out));
      return tmp; 
    };

    // forming all internal 4-index MO integrals
    void form_4index(std::unique_ptr<double[]>& target) const;
    std::unique_ptr<double[]> form_4index() const;
    void form_4index(std::unique_ptr<double[]>& target, const std::shared_ptr<const DF_Full> o) const;
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DF_Full> o) const;

    // the same as above, but with one index fixed at n (of nocc2).
    void form_4index(std::unique_ptr<double[]>& target, const std::shared_ptr<const DF_Full> o, const size_t n) const;
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DF_Full> o, const size_t n) const;

    // forming two-external K operators
    void form_4index(std::unique_ptr<double[]>& target, const std::shared_ptr<const DensityFit> o) const;
    std::unique_ptr<double[]> form_4index(const std::shared_ptr<const DensityFit> o) const;

    double* data() { return data_.get(); };
    const double* const data() const { return data_.get(); };
    const std::shared_ptr<const DensityFit> df() const { return df_; };

};

#endif

