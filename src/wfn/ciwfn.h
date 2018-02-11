//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ciwfn.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef _BAGEL_WFN_CIWFN_H
#define _BAGEL_WFN_CIWFN_H

#include <src/wfn/coeff.h>
#include <src/wfn/geometry.h>
#include <src/ci/fci/dvec.h>
#include <src/ci/zfci/reldvec.h>
#include <src/ci/zfci/relspace.h>

// Stores the result of some CI type wavefunction (FCI, CASSCF, etc.)

namespace bagel {

template<class CIVecClass, class DetClass>
class CIWfn_ {
  protected:
    std::shared_ptr<const Geometry> geom_;

    int ncore_;
    int nact_;

    int nstates_;
    std::shared_ptr<const DetClass> det_;
    std::shared_ptr<const CIVecClass> ccvec_;
    std::vector<double> energies_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & geom_ & ncore_ & nact_ & nstates_ & det_ & ccvec_ & energies_;
    }

  public:
    CIWfn_() { }
    CIWfn_(std::shared_ptr<const Geometry> g, const int ncore, const int nact, const int nst,
          std::vector<double> en, std::shared_ptr<const CIVecClass> ccvec, std::shared_ptr<const DetClass> det)
      : geom_(g), ncore_(ncore), nact_(nact), nstates_(nst),
        det_(det), ccvec_(ccvec), energies_(en) {}

    std::shared_ptr<const Geometry> geom() const { return geom_; }

    int ncore() const { return ncore_; }
    int nact() const { return nact_; }

    int nstates() const { return nstates_; }

    std::vector<double> energies() const { return energies_; }
    double energy(int i) const {return energies_[i];}

    // function to return a CI vectors from orbital info
    std::shared_ptr<const DetClass> det() const { return det_; }
    std::shared_ptr<const CIVecClass> civectors() const { return ccvec_; }
};

using CIWfn = CIWfn_<Dvec,Determinants>;
using RelCIWfn = CIWfn_<RelZDvec,std::pair<std::shared_ptr<const RelSpace>,std::shared_ptr<const RelSpace>>>;

}

#endif
