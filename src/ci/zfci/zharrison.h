//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison.h
// Copyright (C) 2013 Toru Shiozaki
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

// Desc :: The implementation closely follows Harrison and Zarrabian
//

#ifndef __BAGEL_ZFCI_ZHARRISON_H
#define __BAGEL_ZFCI_ZHARRISON_H

#include <src/util/kramers.h>
#include <src/util/math/davidson.h>
#include <src/wfn/method.h>
#include <src/wfn/relreference.h>
#include <src/ci/zfci/zmofile.h>
#include <src/ci/zfci/reldvec.h>

namespace bagel {

class ZHarrison : public Method {

  protected:
    // max #iteration
    int max_iter_;
    int davidson_subspace_;

    // threshold for variants
    double thresh_;
    double print_thresh_;

    // numbers of electrons
    int nele_;
    int ncore_;
    int norb_;
    int charge_;

    // number of states
    int nstate_;
    std::vector<int> states_;

    // total energy
    std::vector<double> energy_;

    // CI vector at convergence
    std::shared_ptr<RelZDvec> cc_;

    // MO integrals
    std::shared_ptr<ZMOFile> jop_;

    // Determinant space
    std::shared_ptr<const RelSpace> space_;
    std::shared_ptr<const RelSpace> int_space_;

    // denominator
    std::shared_ptr<RelDvec> denom_;

    // RDMs
    std::vector<std::shared_ptr<Kramers<2,ZRDM<1>>>> rdm1_;
    std::vector<std::shared_ptr<Kramers<4,ZRDM<2>>>> rdm2_;
    // state averaged RDMs
    std::shared_ptr<Kramers<2,ZRDM<1>>> rdm1_av_;
    std::shared_ptr<Kramers<4,ZRDM<2>>> rdm2_av_;
    std::shared_ptr<ZRDM<1>> rdm1_av_expanded_;
    std::shared_ptr<ZRDM<2>> rdm2_av_expanded_;

    std::shared_ptr<DavidsonDiag<RelZDvec, ZMatrix>> davidson_;

    // integral reuse
    bool store_half_ints_;
    bool store_gaunt_half_ints_;

    // restart
    bool restart_;
    bool restarted_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Method>(*this);
      ar << max_iter_ << davidson_subspace_ << thresh_ << print_thresh_ << nele_ << ncore_ << norb_ << charge_
         << nstate_ << states_ << energy_ << cc_ << space_ << int_space_ << denom_ << rdm1_ << rdm2_ << rdm1_av_ << rdm2_av_ << davidson_
         << store_half_ints_ << store_gaunt_half_ints_ << restart_ << restarted_;
      // for jop_
      std::shared_ptr<const ZCoeff_Block> coeff = jop_->coeff();
      ar << coeff;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Method>(*this);
      ar >> max_iter_ >> davidson_subspace_ >> thresh_ >> print_thresh_ >> nele_ >> ncore_ >> norb_ >> charge_
         >> nstate_ >> states_ >> energy_ >> cc_ >> space_ >> int_space_ >> denom_ >> rdm1_ >> rdm2_ >> rdm1_av_ >> rdm2_av_ >> davidson_
         >> store_half_ints_ >> store_gaunt_half_ints_ >> restart_ >> restarted_;
      std::shared_ptr<const ZCoeff_Block> coeff;
      ar >> coeff;
      update(coeff);
      restarted_ = true;
    }

  protected:
    // obtain determinants for guess generation
    void generate_guess(const int nelea, const int neleb, const int nstate, std::shared_ptr<RelZDvec> out, const int offset);
    // generate spin-adapted guess configurations
    std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet, const int nelea, const int neleb) const;

    // pure virtual functions to be implemented by derived classes
    virtual std::shared_ptr<const ZCoeff_Block> init_coeff() = 0;
    virtual void dump_integrals_and_exit() const = 0;

    void const_denom();

    // run-time functions.
    // aaaa and bbbb
    void sigma_aa(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZMOFile> jop, const bool trans = false) const;
    void sigma_1e_ab(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZMOFile> jop, const bool trans = false) const;

    void sigma_1e_annih_a(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;
    void sigma_1e_annih_b(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;

    void sigma_2e_annih_ab(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;
    void sigma_2e_annih_aa(std::shared_ptr<const ZCivec> cc, std::shared_ptr<ZDvec> d) const;

    void sigma_2e_create_bb(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;
    void sigma_2e_create_ab(std::shared_ptr<ZCivec> sigma, std::shared_ptr<const ZDvec> e) const;

    void sigma_2e_h0101_h1001(std::shared_ptr<const ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const ZMOFile> jop) const;

    template<int i, int j, int k, int l,
             class = typename std::enable_if<i/2 == 0 and j/2 == 0 and k/2 == 0 and l/2 == 0>::type
            >
    void sigma_2e_h(std::shared_ptr<const ZDvec> d, std::shared_ptr<ZDvec> e, std::shared_ptr<const ZMOFile> jop, const bool trans, const std::complex<double> fac = 1.0) const {
      std::bitset<4> bit4((i<<3)+(j<<2)+(k<<1)+l);
      btas::contract(fac, *d, {0,1,2}, *jop->mo2e(trans ? ~bit4 : bit4), {3,2}, std::complex<double>(0.0), *e, {0,1,3});
    }

    void sigma_one(std::shared_ptr<const ZCivec> cc, std::shared_ptr<RelZDvec> sigmavec, std::shared_ptr<const ZMOFile> jop,
                   const int istate, const bool diag, const bool trans) const;
#ifdef HAVE_MPI_H
    int sigma_one_parallel(const int icnt, std::shared_ptr<const ZCivec> cc, std::shared_ptr<RelZDvec> sigmavec, std::shared_ptr<const ZMOFile> jop,
                           const int istate, const bool diag, const bool trans) const;
#endif

    // protected functions for RDM computation. Annihilate two electrons from the reference
    std::shared_ptr<Kramers<1,ZDvec>> one_down_from_civec(const int nelea, const int neleb, const int istate, std::shared_ptr<const RelSpace>) const;
    std::shared_ptr<Kramers<2,ZDvec>> two_down_from_civec(const int nelea, const int neleb, const int istate) const;
    std::shared_ptr<Kramers<3,ZDvec>> three_down_from_civec(const int nelea, const int neleb, const int istate, std::shared_ptr<const RelSpace>) const;
    std::shared_ptr<Kramers<4,ZDvec>> four_down_from_civec(const int nelea, const int neleb, const int istate, std::shared_ptr<const RelSpace>) const;

  public:
    ZHarrison() { }
    // this constructor is ugly... to be fixed some day...
    ZHarrison(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
              const int ncore, const int nocc, std::shared_ptr<const ZCoeff_Block> coeff_zcas, const bool store_c, const bool store_g);

    std::shared_ptr<RelZDvec> form_sigma(std::shared_ptr<const RelZDvec> c, std::shared_ptr<const ZMOFile> jop, const std::vector<int>& conv) const;

    virtual void update(std::shared_ptr<const ZCoeff_Block> coeff) = 0;
    virtual void compute() override;

    // returns members
    int norb() const { return norb_; }
    int ncore() const { return ncore_; }
    int nele() const { return nele_; }
    int nstate() const { return nstate_; }
    double core_energy() const { return jop_->core_energy(); }
    std::shared_ptr<const RelZDvec> cc() const { return cc_; }

    std::shared_ptr<const RelCIWfn> conv_to_ciwfn() const;
    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }

    double energy(const int i) const { return energy_.at(i); }
    std::vector<double> energy() const { return energy_; }

    std::shared_ptr<const ZMOFile> jop() const { return jop_; }

    // functions related to RDMs
    void compute_rdm12();
    std::shared_ptr<Kramers<2,ZRDM<1>>> rdm1(const int jst, const int ist) const;
    std::shared_ptr<Kramers<4,ZRDM<2>>> rdm2(const int jst, const int ist) const;
    std::shared_ptr<Kramers<6,ZRDM<3>>> rdm3(const int jst, const int ist) const;
    std::tuple<std::shared_ptr<Kramers<6,ZRDM<3>>>,std::shared_ptr<Kramers<6,ZRDM<3>>>> rdm34f(const int jst, const int ist, std::shared_ptr<const ZMatrix> fock) const;
    std::shared_ptr<Kramers<8,ZRDM<4>>> rdm4(const int jst, const int ist) const;

    std::vector<std::shared_ptr<const ZMatrix>> rdm1_matrix() const;
    std::shared_ptr<ZMatrix> rdm1_av() const;
    std::shared_ptr<ZMatrix> rdm2_av() const;
    std::shared_ptr<const Kramers<2,ZRDM<1>>> rdm1_av_kramers() const { return rdm1_av_; }
    std::shared_ptr<const Kramers<4,ZRDM<2>>> rdm2_av_kramers() const { return rdm2_av_; }
    template<typename T>
    std::shared_ptr<const ZRDM<1>> rdm1_av_kramers(const T& b) const { KTag<2> t(b); return rdm1_av_->at(t); }
    template<typename T>
    std::shared_ptr<const ZRDM<2>> rdm2_av_kramers(const T& b) const { KTag<4> t(b); return rdm2_av_->at(t); }

    void rotate_rdms(std::shared_ptr<const ZMatrix> trans);

    // interface functions
    void dump_ints() const;
    void read_external_rdm12_av(const std::string& file);
    std::shared_ptr<Kramers<2,ZRDM<1>>> read_external_rdm1(const int ist, const int jst, const std::string& file) const;
    std::shared_ptr<Kramers<4,ZRDM<2>>> read_external_rdm2(const int ist, const int jst, const std::string& file) const;
    std::shared_ptr<Kramers<6,ZRDM<3>>> read_external_rdm3(const int ist, const int jst, const std::string& file, const bool fock_contracted = false) const;
    std::shared_ptr<Kramers<8,ZRDM<4>>> read_external_rdm4(const int ist, const int jst, const std::string& file) const;
};


// only for RDM computation.
class ZFCI_bare : public ZHarrison {
  protected:
    std::shared_ptr<const ZCoeff_Block> init_coeff() override { assert(false); return nullptr; }
    void dump_integrals_and_exit() const override { assert(false); }

  public:
    ZFCI_bare(std::shared_ptr<const RelCIWfn> ci);
    void compute() override { assert(false); }
    void update(std::shared_ptr<const ZCoeff_Block>) override { assert(false); }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ZHarrison)

#endif

