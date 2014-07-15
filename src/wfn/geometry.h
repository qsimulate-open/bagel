//
// BAGEL - Parallel electron correlation program.
// Filename: geometry.h
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


#ifndef __SRC_WFN_GEOMETRY_H
#define __SRC_WFN_GEOMETRY_H

#include <src/df/df.h>
#include <src/input/input.h>
#include <src/molecule/molecule.h>

namespace bagel {

// TODO hide these somewhere
// To enable polymorphic behavior during construction
class Geometry_custom_init {
  public:
    virtual void custom_init() const = 0;
    virtual std::shared_ptr<DFDist> compute_integrals(const double thresh) const = 0;
    virtual int dsize() const = 0;
};

class Geometry_Gaussian_init : public Geometry_custom_init {
  public:
    void custom_init() const override;
    std::shared_ptr<DFDist> compute_integrals(const double thresh) const override;
    int dsize() const override { return 1.0; }
};


class Geometry: public Molecule {
  protected:
    // integral screening
    double schwarz_thresh_;
    double overlap_thresh_;

    // for DF calculations
    mutable std::shared_ptr<DFDist> df_;
    // small component
    mutable std::shared_ptr<DFDist> dfs_;
    // small-large component
    mutable std::shared_ptr<DFDist> dfsl_;

    // Constructor helpers
    void common_init2(const bool print, const double thresh, const bool nodf = false, Geometry_custom_init const& h = Geometry_Gaussian_init());
    void get_electric_field(const std::shared_ptr<const PTree> geominfo);

    // This is to be removed when we eliminate common-origin functionality
    bool london_;
    void set_london(const std::shared_ptr<const PTree> geominfo);

  private:
#if 0
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Molecule>(*this);
      ar << spherical_ << aux_merged_ << nbasis_ << nele_ << nfrc_ << naux_ << lmax_ << aux_lmax_
         << offsets_ << aux_offsets_ << basisfile_ << auxfile_ << schwarz_thresh_ << overlap_thresh_ << london_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Molecule>(*this);
      ar >> spherical_ >> aux_merged_ >> nbasis_ >> nele_ >> nfrc_ >> naux_ >> lmax_ >> aux_lmax_
         >> offsets_ >> aux_offsets_ >> basisfile_ >> auxfile_ >> schwarz_thresh_ >> overlap_thresh_ << london_;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }
#endif

  public:
    Geometry() { }
    Geometry(const std::shared_ptr<const PTree> idata, Geometry_custom_init const& h = Geometry_Gaussian_init());
    Geometry(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o, Geometry_custom_init const& h = Geometry_Gaussian_init());
    Geometry(const Geometry& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true, Geometry_custom_init const& h = Geometry_Gaussian_init());
    Geometry(const Geometry& o, const std::shared_ptr<const Matrix> disp, const std::shared_ptr<const PTree> geominfo, const bool rotate = true, const bool nodf = false);
    Geometry(const Geometry& o, const std::array<double,3> disp);
    Geometry(std::vector<std::shared_ptr<const Geometry>>);

    // Returns a constant
    const std::shared_ptr<const Matrix> compute_grad_vnuc() const;
    double schwarz_thresh() const { return schwarz_thresh_; }
    double overlap_thresh() const { return overlap_thresh_; }
    bool london() const { return london_; }

    // returns schwarz screening TODO not working for DF yet
    std::vector<double> schwarz() const;

    // type T should be a derived class of DFDist
    template<typename T>
    std::shared_ptr<T> form_fit(const double thr, const bool inverse, const double gam = 0.0, const bool average = false) const {
      return std::make_shared<T>(nbasis(), naux(), atoms(), aux_atoms(), thr, inverse, gam, average);
    }

    // Returns DF data
    const std::shared_ptr<const DFDist> df() const { return df_; }
    const std::shared_ptr<const DFDist> dfs() const { return dfs_; }
    const std::shared_ptr<const DFDist> dfsl() const { return dfsl_; }

    // TODO resolve "mutable" issues
    void discard_df() const { df_.reset(); dfs_.reset(); dfsl_.reset(); }

    // initialize relativistic components
    std::shared_ptr<const Geometry> relativistic(const bool do_gaunt) const;
    void compute_relativistic_integrals(const bool do_gaunt);
    void discard_relativistic() const;

};


}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry)

#endif
