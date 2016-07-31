//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_base.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_CI_FCI_FCI_BASE_H
#define __SRC_CI_FCI_FCI_BASE_H

#include <src/ci/fci/mofile.h>
#include <src/wfn/ciwfn.h>
#include <src/wfn/method.h>
#include <src/ci/fci/dvec.h>
#include <src/ci/fci/distcivec.h>
#include <src/util/math/davidson.h>

namespace bagel {

template<class CivecType, class DvecType>
class FCI_base : public Method {
  protected:

    // Options
    int max_iter_;
    int davidson_subspace_;
    int nguess_;
    double thresh_;
    double print_thresh_;

    int nelea_;
    int neleb_;
    int ncore_;
    int norb_;

    int nstate_;

    // extra
    std::shared_ptr<const Determinants> det_;

    // results
    std::vector<double> energy_;
    std::shared_ptr<DvecType> cc_;
    std::shared_ptr<MOFile> jop_;
    // denominator
    std::shared_ptr<CivecType> denom_;

    // RDMs; should be resized in constructors
    std::shared_ptr<VecRDM<1>> rdm1_;
    std::shared_ptr<VecRDM<2>> rdm2_;
    // state averaged RDM
    std::vector<double> weight_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;

    // davidson
    std::shared_ptr<DavidsonDiag<CivecType>> davidson_;

    virtual void print_header() const = 0;

    // const_denom function here only makes a denom for local data of DistCivec.
    virtual void const_denom() = 0;

    // restart
    bool restart_;
    bool restarted_;

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

  public:
    // this constructor is ugly... to be fixed some day...
    FCI_base(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r,
             const int ncore = -1, const int norb = -1, const int nstate = -1)
      : Method(idat, g, r), ncore_(ncore), norb_(norb), nstate_(nstate), restarted_(false) {
    }

    FCI_base() { }
    virtual ~FCI_base() { }

    // FCI compute function
    virtual void compute() override = 0;

    int norb() const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    int ncore() const { return ncore_; }
    double core_energy() const { return jop_->core_energy(); }
    double weight(const int i) const { return weight_.at(i); }

    int nij() const { return (norb_*(norb_+1))/2; }

    virtual void update(std::shared_ptr<const Matrix>) = 0;

    std::shared_ptr<const Determinants> det() const { return det_; }
    std::shared_ptr<const MOFile> jop() const { return jop_; }
    std::shared_ptr<const CivecType> denom() const { return denom_; }
    std::shared_ptr<const DvecType> civectors() const { return cc_; }

    std::vector<double> energy() const { return energy_; }
    double energy(const int i) const { return energy_.at(i); }

    // TODO
    virtual std::pair<std::shared_ptr<Matrix>, VectorB> natorb_convert() = 0;

    virtual std::shared_ptr<const CIWfn> conv_to_ciwfn() const = 0;
    virtual std::shared_ptr<const Reference> conv_to_ref() const = 0; 
};

extern template class FCI_base<Civec,Dvec>;
extern template class FCI_base<DistCivec,DistDvec>;

}

//#include <src/util/archive.h>
//BOOST_CLASS_EXPORT_KEY(bagel::FCI_base<Civec,Dvec>)
//BOOST_CLASS_EXPORT_KEY(bagel::FCI_base<DistCivec,DistDvec>)


#endif
