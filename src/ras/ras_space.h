//
// BAGEL - Parallel electron correlation program.
// Filename: ras_space.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_RAS_RAS_SPACE_H
#define __SRC_RAS_RAS_SPACE_H

#include <src/ras/determinants.h>

namespace bagel {

/// Contains a collection of RASDeterminants to be used in a product space RAS routine
class RASSpace {
  protected:
    std::map<std::pair<int, int>, std::shared_ptr<RASDeterminants>> detmap_; ///< Map to retrieve desired determinants

    std::array<int, 3> ras_;
    int max_holes_;
    int max_particles_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & detmap_;
    }

  public:
    /// Default constructor
    RASSpace() { }
    /** Constructor sets up the basic information on the space
      \param ras number of orbitals in each subspace (RASI, RASII, RASIII)
      \param max_holes maximum number of holes in RASI
      \param max_particles maximum number of particles in RASIII
    */
    RASSpace(const std::array<int, 3> ras, const int max_holes, const int max_particles) : ras_(ras), max_holes_(max_holes), max_particles_(max_particles) { }
    virtual ~RASSpace() { }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    /// Returns a RASDeterminants with the desired number of alpha and beta electrons.
    /// If such a determinant already exists, then it is returned from the map,
    /// if no such determinant already exists, then it is created and returned.
    std::shared_ptr<RASDeterminants> det(const int na, const int nb) {
      auto iter = detmap_.find({na,nb});
      std::shared_ptr<RASDeterminants> out;
      if (iter != detmap_.end()) {
        out = *iter;
      }
      else {
        out = std::make_shared<RASDeterminants>(ras_, na, nb, max_holes_, max_particles_, /*mute*/true);
        detmap_.emplace({na,nb}, out);
      }
      return out;
    }
    /// Const version of det() function only returns a RASDeterminants object if it's already in the detmap_
    std::shared_ptr<const RASDeterminants> det(const int na, const int nb) const { return detmap_.at({na,nb}); }

    const std::map<std::pair<int,int>, std::shared_ptr<RASDeterminants>>& detmap() const { return detmap_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RASSpace)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<RASSpace, T>::value>::type> {
    typedef RASSpace type;
  };
}

#endif
