//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras_space.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __SRC_RAS_RAS_SPACE_H
#define __SRC_RAS_RAS_SPACE_H

#include <src/ci/ras/determinants.h>

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
        out = iter->second;
      }
      else {
        out = std::make_shared<RASDeterminants>(ras_, na, nb, max_holes_, max_particles_, /*mute*/true);
        detmap_.emplace(std::make_pair(na,nb), out);
      }
      return out;
    }
    /// Const version of det() function only returns a RASDeterminants object if it's already in the detmap_
    std::shared_ptr<const RASDeterminants> det(const int na, const int nb) const { return detmap_.at({na,nb}); }

    const std::map<std::pair<int,int>, std::shared_ptr<RASDeterminants>>& detmap() const { return detmap_; }

    bool operator==(const RASSpace& o) { return (ras_==o.ras_ && max_holes_==o.max_holes_ && max_particles_==o.max_particles_); }

    int norb() const { return ras_[0]+ras_[1]+ras_[2]; }
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
