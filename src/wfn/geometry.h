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
#include <src/wfn/geometry_base.h>
#include <src/input/input.h>

namespace bagel {

class Geometry: public Geometry_base {
  protected:

    // for DF calculations
    mutable std::shared_ptr<DFDist> df_;
    // small component
    mutable std::shared_ptr<DFDist> dfs_;
    // small-large component
    mutable std::shared_ptr<DFDist> dfsl_;

    // if false, use common origin with Gaussian orbitals
    bool london_;

    // Constructor helpers
    void compute_integrals(const double thresh, const bool nodf) override;
    void custom_init() override { assert(!nonzero_magnetic_field()); }

    double dsize() const override { return 1.0; }

  private:

    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Geometry_base>(*this);
      const size_t dfindex = !df_ ? 0 : std::hash<DFDist*>()(df_.get());
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
      static std::map<size_t, std::weak_ptr<DFDist>> dfmap;
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
    Geometry() { }
    Geometry(const std::shared_ptr<const PTree>);
    Geometry(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o);
    Geometry(const Geometry& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);
    Geometry(const Geometry& o, const std::shared_ptr<const Matrix> disp, const std::shared_ptr<const PTree> geominfo, const bool rotate = true, const bool nodf = false);
    Geometry(const Geometry& o, const std::array<double,3> disp);
    Geometry(std::vector<std::shared_ptr<const Geometry>>);

    // Returns DF data
    const std::shared_ptr<const DFDist> df() const { return df_; }
    const std::shared_ptr<const DFDist> dfs() const { return dfs_; }
    const std::shared_ptr<const DFDist> dfsl() const { return dfsl_; }

    // TODO resolve "mutable" issues
    void discard_df() const { df_.reset(); dfs_.reset(); dfsl_.reset(); dfints_ = false; }

    // initialize relativistic components
    std::shared_ptr<const Geometry> relativistic(const bool do_gaunt) const;
    void compute_relativistic_integrals(const bool do_gaunt);
    void discard_relativistic() const;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Geometry)

#endif
