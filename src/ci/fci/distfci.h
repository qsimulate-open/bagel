//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distfci.h
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

#ifndef __BAGEL_FCI_DISTFCI_H
#define __BAGEL_FCI_DISTFCI_H

#include <src/ci/fci/distcivec.h>
#include <src/ci/fci/space.h>
#include <src/ci/fci/fci_base.h>

namespace bagel {

// Parallel FCI based on Harrison-Zarrabian algorithm.
// The algorithm is basically the same as written in
// Z. Gan and R. J. Harrison, SC '05: Proceedings of the 2005 ACM/IEEE conference on Supercomputing
//
// The implementation is based on the HarrisonZarrabian class written by Shane Parker.

class DistFCI : public FCI_base {
  protected:
    std::shared_ptr<Space_base> space_;
    std::shared_ptr<DistDvec> cc_;
    std::shared_ptr<DistCivec> denom_;
    std::shared_ptr<DavidsonDiag<DistCivec>> davidson_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Method>(*this);
      ar << max_iter_ << davidson_subspace_ << nguess_ << thresh_ << print_thresh_
         << nelea_ << neleb_ << ncore_ << norb_ << nstate_ << det_
         << energy_ << cc_ << rdm1_ << rdm2_ << weight_ << rdm1_av_ << rdm2_av_ << davidson_;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      // jop_ and denom_ will be constructed in derived classes
      ar >> boost::serialization::base_object<Method>(*this);
      ar >> max_iter_ >> davidson_subspace_ >> nguess_ >> thresh_ >> print_thresh_
         >> nelea_ >> neleb_ >> ncore_ >> norb_ >> nstate_ >> det_
         >> energy_ >> cc_ >> rdm1_ >> rdm2_ >> weight_ >> rdm1_av_ >> rdm2_av_ >> davidson_;
      restarted_ = true;
    }

    void common_init();
    void print_header() const override;

    // const_denom function here only makes a denom for local data of DistCivec.
    void const_denom() override;

    // denominator
    void generate_guess(const int nspin, const int nstate, std::vector<std::shared_ptr<DistCivec>>& out);
    void model_guess(std::vector<std::shared_ptr<DistCivec>>& out);

    // Determinant seeds in parallel
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet) const;

  public:
    // this constructor is ugly... to be fixed some day...
    DistFCI(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
            const int ncore = -1, const int nocc = -1, const int nstate = -1, const bool store = false)
      : FCI_base(a, g, b, ncore, nocc, nstate, store) {
      common_init();
      update(b->coeff());
    }

    // FCI compute function using DistCivec
    void compute() override;
    std::shared_ptr<const Civec> denom() const override;
    std::shared_ptr<const Dvec> civectors() const override;

    void update(std::shared_ptr<const Matrix>) override;

    void compute_rdm12() override;
    void compute_rdm12(const int ist, const int jst) override;
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34(const int ist, const int jst) const override;
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<3>>> rdm34f(const int ist, const int jst, std::shared_ptr<const Matrix> fock) const override;
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>> rdm12_alpha(const int ist, const int jst) const override;
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34_alpha(const int ist, const int jst) const override;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_from_civec(std::shared_ptr<const DistCivec> cbra, std::shared_ptr<const DistCivec> cket) const;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_av_from_dvec(std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>, std::shared_ptr<const Determinants> o) const override;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_last_step(std::shared_ptr<const DistDvec>, std::shared_ptr<const DistDvec>, std::shared_ptr<const DistCivec>) const;

    void sigma_2a1(std::shared_ptr<const DistCivec> c, std::shared_ptr<DistDvec> d) const;
    void sigma_2a2(std::shared_ptr<const DistCivec> c, std::shared_ptr<DistDvec> d) const;

    // rdm ci derivatives
    std::shared_ptr<Dvec> rdm1deriv(const int istate) const override;
    std::shared_ptr<Dvec> rdm2deriv(const int istate) const override;
    std::shared_ptr<Matrix> rdm2fderiv(const int istate, std::shared_ptr<const Matrix> fock, std::shared_ptr<const Matrix> dmat) const override;
    std::shared_ptr<Matrix> rdm2deriv_offset(const int istate, const size_t dsize, const size_t offset, std::shared_ptr<const Matrix> dmat, const bool parallel = true) const override;
    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
      rdm3deriv(const int istate, std::shared_ptr<const Matrix> fock, const size_t offset, const size_t size, std::shared_ptr<const Matrix> dbra_in, std::shared_ptr<const Matrix> fock_ebra_in) const override;

    std::shared_ptr<const CIWfn> conv_to_ciwfn() const override;
    std::shared_ptr<const Dvec> conv_to_dvec() const;
    std::shared_ptr<Dvec> distdvec_to_dvec(std::shared_ptr<const DistDvec> d) const;
    std::shared_ptr<const DistDvec> dvec_to_distdvec(std::shared_ptr<const Dvec> c) const;
    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }

    void read_external_rdm12_av(const std::string& file) override { };
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::DistFCI)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<DistFCI, T>::value>::type> {
    typedef DistFCI type;
  };
}

#endif
