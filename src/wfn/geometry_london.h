//
// BAGEL - Parallel electron correlation program.
// Filename: geometry_london.h
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


#ifndef __SRC_WFN_GEOMETRY_LONDON_H
#define __SRC_WFN_GEOMETRY_LONDON_H

#include <src/df/complexdf.h>
#include <src/wfn/geometry_base.h>
#include <src/input/input.h>

namespace bagel {

class Geometry_London: public Geometry_base {
  protected:

    // for DF calculations
    mutable std::shared_ptr<ComplexDFDist> df_;
    // small component
    mutable std::shared_ptr<ComplexDFDist> dfs_;
    // small-large component
    mutable std::shared_ptr<ComplexDFDist> dfsl_;

    // if false, use common origin with Gaussian orbitals
    bool london_;

    // Constructor helpers
    void compute_integrals(const double thresh, const bool nodf) override;
    void custom_init() override {
      if (london_ && nonzero_magnetic_field()) std::cout << "  Using London orbital basis to enforce gauge-invariance" << std::endl;
      if (!london_ && nonzero_magnetic_field()) std::cout << "  Using a common gauge origin - NOT RECOMMENDED for accurate calculations.  (Use a London orbital basis instead.)" << std::endl;
      if (!nonzero_magnetic_field()) std::cout << "  Zero magnetic field - This computation would be more efficient with a Gaussian basis set." << std::endl;
      magnetic();
    }

    double dsize() const override { return 2.0; }

  private:

    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Geometry_base>(*this);
      const size_t dfindex = !df_ ? 0 : std::hash<ComplexDFDist*>()(df_.get());
      ar << dfindex;
      const bool do_rel   = !!dfs_;
      const bool do_gaunt = !!dfsl_;
      ar << do_rel << do_gaunt;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Geometry_base>(*this);
      size_t dfindex;
      ar >> dfindex;
      static std::map<size_t, std::weak_ptr<ComplexDFDist>> dfmap;
      if (dfmap[dfindex].expired()) {
        compute_integrals(overlap_thresh_, dfindex == 0);
        dfmap[dfindex] = df_;
      } else {
        df_ = dfmap[dfindex].lock();
      }

      bool do_rel, do_gaunt;
      ar >> do_rel >> do_gaunt;
      if (do_rel)
        compute_relativistic_integrals(do_gaunt);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    Geometry_London() { }
    Geometry_London(const std::shared_ptr<const PTree>);
    Geometry_London(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o);
    Geometry_London(const Geometry_London& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);

    // Returns a constant
    bool london() const {return london_; }

    // Returns DF data
    const std::shared_ptr<const ComplexDFDist> df() const { return df_; }
    const std::shared_ptr<const ComplexDFDist> dfs() const { return dfs_; }
    const std::shared_ptr<const ComplexDFDist> dfsl() const { return dfsl_; }

    // TODO resolve "mutable" issues
    void discard_df() const { df_.reset(); dfs_.reset(); dfsl_.reset(); }

    // initialize relativistic components
    std::shared_ptr<const Geometry_London> relativistic(const bool do_gaunt) const;
    void compute_relativistic_integrals(const bool do_gaunt);
    void discard_relativistic() const;

    // initialize magnetic components
    void magnetic();

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry_London)

#endif
