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
#include <src/wfn/geometry.h>
#include <src/input/input.h>

namespace bagel {

class Geometry_London_init : public Geometry_custom_init {
    void custom_init() const override;
    std::shared_ptr<DFDist> compute_integrals(const double thresh) const override;
    int dsize() const override { return 2.0; }
};


class Geometry_London: public Geometry {
  protected:

  private:

#if 0
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Geometry>(*this);
      const size_t dfindex = !df_ ? 0 : std::hash<ComplexDFDist*>()(df_.get());
      ar << dfindex;
      const bool do_rel   = !!dfs_;
      const bool do_gaunt = !!dfsl_;
      ar << do_rel << do_gaunt;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Geometry>(*this);
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
#endif

  public:
    Geometry_London() { }
    Geometry_London(const std::shared_ptr<const PTree> idata);
    Geometry_London(const std::vector<std::shared_ptr<const Atom>> atoms, const std::shared_ptr<const PTree> o);
    Geometry_London(const Geometry_London& o, const std::shared_ptr<const PTree> idata, const bool discard_prev_df = true);

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
