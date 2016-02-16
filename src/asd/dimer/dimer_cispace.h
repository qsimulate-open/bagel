//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer_cispace.h
// Copyright (C) 2012 Toru Shiozaki
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



#ifndef __BAGEL_DIMER_CISPACE_H
#define __BAGEL_DIMER_CISPACE_H

#include <utility>
#include <set>
#include <src/ci/fci/dvec.h>
#include <src/ci/ras/civector.h>

namespace bagel {

class SpaceKey {
  public:
    const int S;
    const int m_s;
    const int q;

    SpaceKey(const int _s, const int _m_s, const int _q) : S(_s), m_s(_m_s), q(_q) {}

    bool operator<(const SpaceKey& o) const {
      if ( abs(q) != abs(o.q)) return ( abs(q) < abs(o.q));
      else if ( q != o.q) return (q < o.q);
      else if (S != o.S) return (S < o.S);
      else return (m_s < o.m_s);
    }
    bool operator==(const SpaceKey& o) const { return ( ( (S == o.S) && (m_s == o.m_s) ) && (q == o.q) ); }

    std::string to_string() const {
      std::vector<std::string> anions = {{ "O", "A", "diA", "triA", "tetA", "pentA" }};
      std::vector<std::string> cations = {{ "O", "C", "diC", "triC", "tetC", "pentC" }};

      std::string out = ( q < 0 ? anions.at(-q) : cations.at(q) );
      out = out + "(2S=" + std::to_string(S) + ";2m_s=" + std::to_string(m_s) + ")";
      return out;
    }

    // unique tag
    int tag() const { return (((S << 5) + m_s+S) << 5) + q; }
};

template <class VecType>
class DimerCISpace_base {
  using DetType = typename VecType::DetType;

  using SpaceMap = std::map<SpaceKey, std::shared_ptr<VecType>>;
  using DMap = std::map<std::pair<int,int>, std::shared_ptr<DetType>>;

  protected:
    // These are stored values of the neutral species
    std::pair<int, int> nelea_;
    std::pair<int, int> neleb_;
    std::pair<int, int> nstates_;

    // This is used only to make new DetType objects through clone()
    std::pair<std::shared_ptr<const DetType>, std::shared_ptr<const DetType>> bdet_;

    SpaceMap cispaceA_;
    SpaceMap cispaceB_;

    DMap detspaceA_;
    DMap detspaceB_;


  public:
    // This constructor will build the infrastructure; civecs need to be added later
    DimerCISpace_base(std::pair<std::shared_ptr<const DetType>, std::shared_ptr<const DetType>> bdet, std::pair<int, int> nelea, std::pair<int, int> neleb) : nelea_(nelea), neleb_(neleb), bdet_(bdet) {}

    template<int unit> int nelea() const { return (unit == 0 ? nelea_.first : nelea_.second); }
    template<int unit> int neleb() const { return (unit == 0 ? neleb_.first : neleb_.second); }
    template<int unit> int nstates() const { return (unit == 0 ? nstates_.first : nstates_.second); }
    template<int unit> std::shared_ptr<const DetType> bdet() const { return (unit == 0 ? bdet_.first : bdet_.second); }
    template<int unit> int norb() const { return unit == 0 ? detspaceA_.cbegin()->second->norb() : detspaceB_.cbegin()->second->norb(); }

    std::pair<int, int> nelea() const { return nelea_; }
    std::pair<int, int> neleb() const { return neleb_; }
    std::pair<int, int> nstates() const { return nstates_; }

    template<int unit> SpaceMap& cispace() { return (unit == 0 ? cispaceA_ : cispaceB_); }
    template<int unit> const SpaceMap& cispace() const { return (unit == 0 ? cispaceA_ : cispaceB_); }

    template<int unit> std::shared_ptr<VecType> ccvec(const int S, const int m_s, const int q) { return ccvec<unit>(SpaceKey(S,m_s,q)); }
    template<int unit> std::shared_ptr<VecType> ccvec(SpaceKey key) {
      SpaceMap& space = (unit == 0 ? cispaceA_ : cispaceB_);
      auto iter = space.find(key);
      return (iter != space.end() ? iter->second : nullptr);
    }

    template<int unit> std::set<SpaceKey> spacekeys() const {
      std::set<SpaceKey> out;
      for (auto& i : cispace<unit>()) out.insert(i.first);
      return out;
    }

    template<int unit> std::shared_ptr<DetType> det(std::pair<const int, const int> p) { return det<unit>(p.first,p.second); }
    template<int unit> std::shared_ptr<DetType> det(const int qa, const int qb) {
      DMap& dets = (unit == 0 ? detspaceA_ : detspaceB_);
      auto iter = dets.find({qa, qb});
      return (iter != dets.end() ? iter->second : nullptr);
    }

    template<int unit, class T> void insert(std::shared_ptr<const T> civec, const int spin = -1) {
      auto new_civec = std::make_shared<VecType>(*civec);

      const int nelea = civec->det()->nelea();
      const int neleb = civec->det()->neleb();

      const int m_s = nelea - neleb;
      const int S = (spin < 0) ? m_s : spin;
      const int Q = charge<unit>(nelea, neleb);

      // Reform DetType object (to make sure it's the format I want)
      std::shared_ptr<DetType> det = add_det<unit>(detkey<unit>(nelea, neleb));
      new_civec->set_det(det);

      SpaceMap& cispace = (unit == 0 ? cispaceA_ : cispaceB_);
      cispace.emplace(SpaceKey(S,m_s,Q), new_civec);

      const int ij = new_civec->ij();
      ((unit == 0) ? nstates_.first : nstates_.second) += ij;
    }

    template<int unit> std::shared_ptr<DetType> add_det(const int qa, const int qb) {
      DMap& detspace = (unit == 0 ? detspaceA_ : detspaceB_);

      typename DMap::iterator idet = detspace.find({qa,qb});
      if ( idet != detspace.end()) {
        return idet->second;
      }
      else {
        int nelea, neleb;
        std::tie(nelea, neleb) = detunkey<unit>(qa,qb);
        std::shared_ptr<DetType> det = bdet<unit>()->clone(nelea, neleb);

        detspace.emplace(std::make_pair(qa,qb), det);

        idet = detspace.find({qa+1,qb});
        if (idet != detspace.end()) det->template link<0>(idet->second);

        idet = detspace.find({qa-1,qb});
        if (idet != detspace.end()) det->template link<0>(idet->second);

        idet = detspace.find({qa,qb+1});
        if (idet != detspace.end()) det->template link<1>(idet->second);

        idet = detspace.find({qa,qb-1});
        if (idet != detspace.end()) det->template link<1>(idet->second);

        return det;
      }
    }

    template<int unit> std::shared_ptr<DetType> add_det(std::pair<const int, const int> p) { return add_det<unit>(p.first,p.second); }

    void insert(std::pair<std::shared_ptr<const VecType>, std::shared_ptr<const VecType>> cipair) { insert<0>(cipair.first); insert<1>(cipair.second); }

    // Completes the spin case and adds extra Determinants into the mapping that will be needed for Hamiltonian computation
    void complete() {
      {
        std::vector<SpaceKey> references;
        for (auto& imap : cispaceA_) { references.push_back(imap.first); }
        // These spaces are assumed to be high-spin
        for (auto& ispace : references) {
          if (ispace.S > 0) {
            const int S = ispace.S;

            std::shared_ptr<VecType> ref_state = ccvec<0>(ispace);

            int ref_qa, ref_qb;
            std::tie(ref_qa, ref_qb) = detkey<0>(ispace);
            const int mult = S + 1;

            for (int i = 1; i < mult; ++i) {
              const int nqa = ref_qa + i;
              const int nqb = ref_qb - i;

              std::shared_ptr<DetType> det = add_det<0>(nqa, nqb);

              ref_state = ref_state->spin_lower(det);
              for (int istate = 0; istate < ref_state->ij(); ++istate) {
                const double norm = ref_state->data(istate)->norm();
                if ( norm < numerical_zero__ ) throw std::runtime_error("Spin lowering operator yielded no state.");
                ref_state->data(istate)->scale(1.0/norm);
              }
              insert<0>(std::shared_ptr<const VecType>(ref_state), S);
            }
          }
        }
      }


      {
        std::vector<SpaceKey> references;
        for (auto& imap : cispaceB_) { references.push_back(imap.first); }
        // These spaces are assumed to be high-spin
        for (auto& ispace : references) {
          if (ispace.S > 0) {
            const int S = ispace.S;

            std::shared_ptr<VecType> ref_state = ccvec<1>(ispace);

            int ref_qa, ref_qb;
            std::tie(ref_qa, ref_qb) = detkey<1>(ispace);
            const int mult = S + 1;

            for (int i = 1; i < mult; ++i) {
              const int nqa = ref_qa + i;
              const int nqb = ref_qb - i;

              std::shared_ptr<DetType> det = add_det<1>(nqa, nqb);

              ref_state = ref_state->spin_lower(det);
              for (int istate = 0; istate < ref_state->ij(); ++istate) {
                const double norm = ref_state->data(istate)->norm();
                if ( norm < numerical_zero__ ) throw std::runtime_error("Spin lowering operator yielded no state.");
                ref_state->data(istate)->scale(1.0/norm);
              }
              insert<1>(std::shared_ptr<const VecType>(ref_state), S);
            }
          }
        }
      }
    }

    // Fills in N-1 electron intermediates for Harrison--Zarrabian CAS algorithm
    void intermediates() {
      for (auto& imap : cispaceA_) {
        int qa, qb;
        std::tie(qa, qb) = detkey<0>(imap.first);

        add_det<0>(qa+1, qb);
        add_det<0>(qa, qb+1);
        add_det<0>(qa+1, qb+1);
      }

      for (auto& imap : cispaceB_) {
        int qa, qb;
        std::tie(qa, qb) = detkey<1>(imap.first);

        add_det<1>(qa+1, qb);
        add_det<1>(qa, qb+1);
        add_det<1>(qa+1, qb+1);
      }
    }

  private:
    template<int unit> std::pair<int, int> detkey(const SpaceKey key) { return detkey<unit>(key.S, key.m_s, key.q); }
    template<int unit> std::pair<int, int> detkey(const int S, const int m_s, const int q) {
      const int nS = m_s - (nelea<unit>() - neleb<unit>());

      return {(q - nS)/2, (q + nS)/2};
    }
    template<int unit> std::pair<int, int> detkey(const int nea, const int neb) const {
      return {nelea<unit>() - nea, neleb<unit>() - neb};
    }

    template<int unit> std::pair<int, int> detunkey(const int qa, const int qb) const {
      return {nelea<unit>() - qa, neleb<unit>() - qb};
    }

    template<int unit> int charge(const int nea, const int neb) const { return ( (nelea<unit>() + neleb<unit>()) - (nea + neb) ); }
};

using DimerCAS = DimerCISpace_base<CASDvec>;
using DimerRAS = DimerCISpace_base<RASDvec>;

}

#endif
