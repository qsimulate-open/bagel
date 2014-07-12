//
// BAGEL - Parallel electron correlation program.
// Filename: geometry_base.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_WFN_GEOMETRY_BASE_H
#define __SRC_WFN_GEOMETRY_BASE_H

#include <src/molecule/molecule.h>

namespace bagel {

class Geometry_base : public Molecule {
  protected:
    // integral screening
    double schwarz_thresh_;
    double overlap_thresh_;

    // Have integrals been calculated?
    bool dfints_;

    // Constructor helpers
    void common_init2(const bool print, const double thresh, const bool nodf = false);
    void get_electric_field(const std::shared_ptr<const PTree> geominfo);

    virtual void compute_integrals(const double thresh, const bool nodf) = 0;
    virtual void custom_init() = 0;

    // So storage requirements are accurately calculated
    virtual double dsize() const = 0;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Molecule>(*this);
      ar << spherical_ << aux_merged_ << nbasis_ << nele_ << nfrc_ << naux_ << lmax_ << aux_lmax_
         << offsets_ << aux_offsets_ << basisfile_ << auxfile_ << schwarz_thresh_ << overlap_thresh_ << dfints_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Molecule>(*this);
      ar >> spherical_ >> aux_merged_ >> nbasis_ >> nele_ >> nfrc_ >> naux_ >> lmax_ >> aux_lmax_
         >> offsets_ >> aux_offsets_ >> basisfile_ >> auxfile_ >> schwarz_thresh_ >> overlap_thresh_ >> dfints_;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    Geometry_base() { }
    Geometry_base(const std::shared_ptr<const PTree>);
    Geometry_base(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o);
    Geometry_base(const Geometry_base& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);
    Geometry_base(const Geometry_base& o, const std::shared_ptr<const Matrix> disp, const std::shared_ptr<const PTree> geominfo, const bool rotate = true, const bool nodf = false);
    Geometry_base(const Geometry_base& o, const std::array<double,3> disp);

    // Returns a constant
    const std::shared_ptr<const Matrix> compute_grad_vnuc() const;
    double schwarz_thresh() const { return schwarz_thresh_; }
    double overlap_thresh() const { return overlap_thresh_; }

    // careful, discarding integrals does not reset dfints_
    double dfints() const { return dfints_; }

    // returns schwarz screening TODO not working for DF yet
    std::vector<double> schwarz() const;

    // type T should be a derived class of DFDist
    template<typename T>
    std::shared_ptr<T> form_fit(const double thr, const bool inverse, const double gam = 0.0, const bool average = false) const {
      return std::make_shared<T>(nbasis(), naux(), atoms(), aux_atoms(), thr, inverse, gam, average);
    }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry_base)

#endif
