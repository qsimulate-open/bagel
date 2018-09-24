//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimersubspace.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_ASD_DIMESUBSPACE_H
#define __SRC_ASD_DIMESUBSPACE_H

#include <src/util/string_util.h>
#include <src/util/math/csymmatrix.h>
#include <src/asd/dimer/dimer_cispace.h>

namespace bagel {

// Space Key + Nstate
class MonomerKey : public SpaceKey {
  protected:
    int nstates_;
    int nelea_;
    int neleb_;
  public:
    MonomerKey(const SpaceKey& s, const int n, const int na, const int nb) : SpaceKey(s), nstates_(n), nelea_(na), neleb_(nb) { }
    int nstates() const { return nstates_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
};

class DimerSubspace_base {
  protected:
    const int offset_;
    std::pair<MonomerKey, MonomerKey> key_;

    std::pair<std::shared_ptr<const CSymMatrix>, std::shared_ptr<const CSymMatrix>> sigma_;

  public:
    DimerSubspace_base(const int o, const MonomerKey& a, const MonomerKey& b) : offset_(o), key_({a,b}) { }

    int nstatesA() const { return key_.first.nstates(); }
    int nstatesB() const { return key_.second.nstates(); }
    template <int unit> int nstates() const { return unit == 0 ? nstatesA() : nstatesB(); }

    std::string stringA() const { return key_.first.to_string(); }
    std::string stringB() const { return key_.second.to_string(); }

    int tagA() const { return key_.first.tag(); }
    int tagB() const { return key_.second.tag(); }
    template <int unit> int tag() const { return unit == 0 ? tagA() : tagB(); }

    std::pair<int,int> S() const { return {key_.first.S, key_.second.S}; }
    std::pair<int,int> ms() const { return {key_.first.m_s, key_.second.m_s}; }
    std::pair<int,int> charge() const { return {key_.first.q, key_.second.q}; }

    int dimerstates() const { return nstatesA() * nstatesB(); }
    int dimerindex(const int iA, const int iB) const { return (iA + iB*nstatesA()); }
    std::string string(const int i, const int j) const {
      std::string out = stringA() + lexical_cast<std::string>(i) + std::string(" ") + stringB() + lexical_cast<std::string>(j);
      return out;
    }

    int offset() const { return offset_; }

    template <int unit> MonomerKey monomerkey() const { return unit == 0 ? MonomerKey(key_.first) : MonomerKey(key_.second); }

    template <int unit> std::shared_ptr<const CSymMatrix> sigma() const { return ( unit == 0 ? sigma_.first : sigma_.second ); }

    template <int unit> void set_sigma(std::shared_ptr<const CSymMatrix> s) { (unit == 0 ? sigma_.first : sigma_.second) = s; }
};


/// Contains all of the information for a product of two monomer spaces
template <class VecType>
class DimerSubspace : public DimerSubspace_base {
  protected:
    std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> ci_;

  public:
    DimerSubspace(int& _offset, const SpaceKey Akey, const SpaceKey Bkey, std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> _ci) :
      DimerSubspace_base(_offset, MonomerKey(Akey, _ci.first->ij(), _ci.first->det()->nelea(), _ci.first->det()->neleb()),
                                  MonomerKey(Bkey, _ci.first->ij(), _ci.second->det()->nelea(), _ci.second->det()->neleb())), ci_({_ci.first, _ci.second})
    { _offset += dimerstates(); }

    template <int unit> std::shared_ptr<const VecType> ci() const { return unit == 0 ? ci_.first : ci_.second; }

};

}

#endif
