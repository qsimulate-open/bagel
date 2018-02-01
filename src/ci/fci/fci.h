//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci.h
// Copyright (C) 2011 Toru Shiozaki
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

#ifndef __BAGEL_FCI_FCI_H
#define __BAGEL_FCI_FCI_H

#include <src/ci/fci/dvec.h>
#include <src/ci/fci/fci_base.h>

namespace bagel {

class FCI : public FCI_base {
  protected:
    std::shared_ptr<Dvec> cc_;
    std::shared_ptr<Civec> denom_;
    std::shared_ptr<DavidsonDiag<Civec>> davidson_;

    bool dipoles_;

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


  protected:
    // some init functions
    void common_init();
    void create_Jiiii();
    // obtain determinants for guess generation
    void generate_guess(const int nspin, const int nstate, std::shared_ptr<Dvec>);
    void model_guess(std::shared_ptr<Dvec>);
    // generate spin-adapted guess configurations
    virtual std::vector<std::pair<std::bitset<nbit__>, std::bitset<nbit__>>> detseeds(const int ndet) const;

    /* Virtual functions -- these MUST be defined in the derived class*/
    // denominator
    virtual void const_denom() override = 0;

    // print functions
    void print_header() const override;

  public:
    FCI() { }

    // this constructor is ugly... to be fixed some day...
    FCI(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>,
        const int ncore = -1, const int nocc = -1, const int nstate = -1, const bool store = false);

    virtual ~FCI() { }
    virtual void compute() override;
    std::shared_ptr<const Civec> denom() const override { return denom_; }
    std::shared_ptr<const Dvec> civectors() const override { return cc_; }

    virtual void update(std::shared_ptr<const Matrix>) = 0;

    void compute_rdm12() override;
    void compute_rdm12(const int ist, const int jst) override;

    // virtual application of Hamiltonian
    virtual std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> c, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const = 0;

    // compute 3 and 4 RDMs
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34(const int ist, const int jst) const override;
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<3>>> rdm34f(const int ist, const int jst, std::shared_ptr<const Matrix> fock) const override;
    // compute "alpha" 1 and 2 RDMs <ia ja> and <ia ja, k, l>
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>> rdm12_alpha(const int ist, const int jst) const override;
    // compute "alpha" 3 and 4 RDMs <ia ja, k, l, m n>...
    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> rdm34_alpha(const int ist, const int jst) const override;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_from_civec(std::shared_ptr<const Civec> cbra, std::shared_ptr<const Civec> cket) const;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_av_from_dvec(std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>, std::shared_ptr<const Determinants> o) const override;

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
      compute_rdm12_last_step(std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>, std::shared_ptr<const Civec>) const;

    // rdm ci derivatives
    std::shared_ptr<Dvec> rdm1deriv(const int istate) const override;
    std::shared_ptr<Dvec> rdm2deriv(const int istate) const override;
    std::shared_ptr<Matrix> rdm2fderiv(const int istate, std::shared_ptr<const Matrix> fock, std::shared_ptr<const Matrix> dmat) const override;
    std::shared_ptr<Matrix> rdm2deriv_offset(const int istate, const size_t dsize, const size_t offset, std::shared_ptr<const Matrix> dmat, const bool parallel = true) const override;
    std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>>
      rdm3deriv(const int istate, std::shared_ptr<const Matrix> fock, const size_t offset, const size_t size, std::shared_ptr<const Matrix> dbra_in, std::shared_ptr<const Matrix> fock_ebra_in) const override;

    // functions for RDM computation
    void sigma_2a1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;
    void sigma_2a2(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;

    std::shared_ptr<const CIWfn> conv_to_ciwfn() const override;
    std::shared_ptr<const Reference> conv_to_ref() const override { return nullptr; }

    // interface functions
    // read state averaged RDM1 and 2 and set to rdm1_av_expanded_ and rdm2_av_expanded_
    void read_external_rdm12_av(const std::string& file) override;
    std::shared_ptr<RDM<1>> read_external_rdm1(const int ist, const int jst, const std::string& file) const;
    std::shared_ptr<RDM<2>> read_external_rdm2(const int ist, const int jst, const std::string& file) const;
    std::shared_ptr<RDM<3>> read_external_rdm3(const int ist, const int jst, const std::string& file, const bool fock_contracted = false) const;
    std::shared_ptr<RDM<4>> read_external_rdm4(const int ist, const int jst, const std::string& file) const;
    void dump_ints() const;
};


// only for RDM computation.
class FCI_bare : public FCI {
  protected:
    void const_denom() override { assert(false); }

  public:
    FCI_bare(std::shared_ptr<const CIWfn> ci);

    void compute() override { assert(false); }
    void update(std::shared_ptr<const Matrix>) override { assert(false); }
    std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec>, std::shared_ptr<const MOFile>, const std::vector<int>&) const override { assert(false); return nullptr; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::FCI)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<FCI, T>::value>::type> {
    typedef FCI type;
  };
}

#endif
