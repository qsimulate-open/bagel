//
// BAGEL - Parallel electron correlation program.
// Filename: dimersubspace.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_MEH_DIMESUBSPACE_H
#define __SRC_MEH_DIMESUBSPACE_H

#include <src/dimer/dimer_cispace.h>

namespace bagel {

// Space Key + Nstate
class MonomerKey : public SpaceKey {
  protected:
    int nstates_;
    int nelea_;
    int neleb_;
  public:
    MonomerKey(const int S, const int ms, const int charge, const int n, const int na, const int nb) : SpaceKey(S,ms,charge), nstates_(n), nelea_(na), neleb_(nb) { }
    MonomerKey(const SpaceKey& s, const int n, const int na, const int nb) : SpaceKey(s), nstates_(n), nelea_(na), neleb_(nb) { }
    int nstates() const { return nstates_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
};

/// Wrapper for monomer CI wavefunctions that includes extra helpful information
template <class VecType>
class MonomerSubspace_base : public MonomerKey {
  protected:
    std::shared_ptr<const VecType> monomerci_;

  public:
    MonomerSubspace_base(const SpaceKey& s, const int nst, std::shared_ptr<const VecType> monomerci) :
      MonomerKey(s, nst, monomerci->det()->nelea(), monomerci->det()->neleb()), monomerci_(monomerci) {}
    MonomerSubspace_base(const int S, const int ms, const int charge, const int nst, std::shared_ptr<const VecType> monomerci) :
      MonomerKey(S, ms, charge, nst,  monomerci->det()->nelea(), monomerci->det()->neleb()), monomerci_(monomerci) {}

    std::shared_ptr<const VecType> monomerci() const { return monomerci_; }
};


/// Contains all of the information for a product of two monomer spaces
template <class VecType>
class DimerSubspace_base {
  protected:
    const int offset_;
    const int nstatesA_;
    const int nstatesB_;
    const std::string stringA_;
    const std::string stringB_;

    std::pair<MonomerSubspace_base<VecType>, MonomerSubspace_base<VecType>> ci_;
    std::pair<std::shared_ptr<const CSymMatrix>, std::shared_ptr<const CSymMatrix>> sigma_;

  public:
    DimerSubspace_base(int& _offset, const SpaceKey Akey, const SpaceKey Bkey, std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> _ci) :
      offset_(_offset), nstatesA_(_ci.first->ij()), nstatesB_(_ci.second->ij()), stringA_(Akey.to_string()), stringB_(Bkey.to_string()),
       ci_({MonomerSubspace_base<VecType>(Akey, nstatesA_, _ci.first), MonomerSubspace_base<VecType>(Bkey, nstatesB_, _ci.second)})
    { _offset += dimerstates(); }

    int offset() const { return offset_; }

    int dimerstates() const { return nstatesA_ * nstatesB_; }
    int dimerindex(const int iA, const int iB) const { return (iA + iB*nstatesA_); }
    std::string string(const int i, const int j) const {
      std::string out = stringA_ + lexical_cast<std::string>(i) + std::string(" ") + stringB_ + lexical_cast<std::string>(j);
      return out;
    }

    template <int unit> int nstates() const { return ( unit == 0 ? nstatesA_ : nstatesB_ ); }
    template <int unit> std::shared_ptr<const VecType> ci() const { return ( unit == 0 ? ci_.first.monomerci() : ci_.second.monomerci() ); }
    template <int unit> std::shared_ptr<const CSymMatrix> sigma() const { return ( unit == 0 ? sigma_.first : sigma_.second ); }
    template <int unit> int tag() const { return unit == 0 ? ci_.first.tag() : ci_.second.tag(); }

    std::pair<int,int> S() const { return {ci_.first.S, ci_.second.S}; }
    std::pair<int,int> ms() const { return {ci_.first.m_s, ci_.second.m_s}; }
    std::pair<int,int> charge() const { return {ci_.first.q, ci_.second.q}; }

    template <int unit> MonomerKey monomerkey() const { return unit == 0 ? MonomerKey(ci_.first) : MonomerKey(ci_.second); }

    template <int unit> void set_sigma(std::shared_ptr<const CSymMatrix> s) { (unit == 0 ? sigma_.first : sigma_.second) = s; }

};

}

#endif
