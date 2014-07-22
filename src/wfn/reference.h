//
// BAGEL - Parallel electron correlation program.
// Filename: reference.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef _BAGEL_WFN_REFERENCE_H
#define _BAGEL_WFN_REFERENCE_H

#include <set>
#include <src/scf/coeff.h>
#include <src/molecule/hcore.h>
#include <src/wfn/geometry.h>
#include <src/fci/dvec.h>
#include <src/wfn/ciwfn.h>
#include <src/wfn/rdm.h>

// all the info to construct wave functions

namespace bagel {

class Reference : public std::enable_shared_from_this<Reference> {

  protected:
    // Geometry which this wave function is belonging to
    std::shared_ptr<const Geometry> geom_;
    // MO coefficients
    std::shared_ptr<const Coeff> coeff_;
    // in case of spin-broken wave functions (UHF)
    std::shared_ptr<const Coeff> coeffA_;
    std::shared_ptr<const Coeff> coeffB_;
    int noccA_, noccB_;

    double energy_;

    std::shared_ptr<const Hcore> hcore_;
    VectorB eig_;

    int nclosed_;
    int nact_;
    int nvirt_;

    int nstate_;

    std::shared_ptr<const CIWfn> ciwfn_;
    std::vector<std::shared_ptr<RDM<1>>>  rdm1_;
    std::vector<std::shared_ptr<RDM<2>>>  rdm2_;
    std::shared_ptr<const RDM<1>> rdm1_av_;
    std::shared_ptr<const RDM<2>> rdm2_av_;

    // this is only for UHF gradient. Somehow I cannot come up with a beautiful design for this.
    std::shared_ptr<const Matrix> erdm1_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & geom_ & coeff_ & coeffA_ & coeffB_ & noccA_ & noccB_ & energy_ & hcore_ & eig_
         & nclosed_ & nact_ & nvirt_ & nstate_ & ciwfn_ & rdm1_ & rdm2_ & rdm1_av_ & rdm2_av_ & erdm1_;
    }

  public:
    Reference() { }
    Reference(std::shared_ptr<const Geometry> g, std::shared_ptr<const Coeff> c,
              const int nclo, const int nact, const int nvirt,
              const double en = 0.0,
              std::vector<std::shared_ptr<RDM<1>>> rdm1 = std::vector<std::shared_ptr<RDM<1>>>(),
              std::vector<std::shared_ptr<RDM<2>>> rdm2 = std::vector<std::shared_ptr<RDM<2>>>(),
              std::shared_ptr<const RDM<1>> rdm1_av = nullptr,
              std::shared_ptr<const RDM<2>> rdm2_av = nullptr,
              std::shared_ptr<const CIWfn> ci = nullptr);

    // copy construct with optionally updating coeff
    Reference(const Reference& o, std::shared_ptr<const Coeff> c = nullptr) :
      Reference(o.geom(), c ? c : o.coeff(), o.nclosed(), o.nact(), o.nvirt(), o.energy(), o.rdm1(), o.rdm2(), o.rdm1_av(), o.rdm2_av(), o.ciwfn()) { }

    virtual ~Reference() { }

    const std::shared_ptr<const Geometry> geom() const { return geom_; }
    const std::vector<double> schwarz() const { return geom_->schwarz(); }
    virtual const std::shared_ptr<const Hcore> hcore() const { return hcore_; }
    virtual const std::shared_ptr<const Coeff> coeff() const { return coeff_; }

    void set_coeff(std::shared_ptr<const Matrix> matrix) { coeff_ = std::make_shared<const Coeff>(*matrix); }
    void set_nocc(const int a, const int b) { noccA_ = a; noccB_ = b; }

    void set_eig(const VectorB& eig);
    const VectorB& eig() const { return eig_; }
    void set_erdm1(const std::shared_ptr<const Matrix> o);
    std::shared_ptr<const Matrix> erdm1() const { return erdm1_; }

    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nvirt() const { return nvirt_; }
    int nocc() const { return nclosed_ + nact_; }

    int noccA() const { return noccA_; }
    int noccB() const { return noccB_; }

    std::shared_ptr<Reference> set_active(std::set<int> active_indices) const;
    std::shared_ptr<Reference> set_ractive(std::set<int> ras1, std::set<int> ras2, std::set<int> ras3) const;

    // used in SA-CASSCF
    void set_nstate(const int i) { nstate_ = i; }
    int nstate() const { return nstate_; }

    // used in UHF
    void set_coeff_AB(const std::shared_ptr<const Coeff> a, const std::shared_ptr<const Coeff> b);
    const std::shared_ptr<const Coeff> coeffA() const { return coeffA_; }
    const std::shared_ptr<const Coeff> coeffB() const { return coeffB_; }

    double energy() const { return energy_; }

    std::shared_ptr<const CIWfn> ciwfn() const { return ciwfn_; }

    const std::vector<std::shared_ptr<RDM<1>>>& rdm1() const { return rdm1_; }
    const std::vector<std::shared_ptr<RDM<2>>>& rdm2() const { return rdm2_; }

    std::shared_ptr<const RDM<1>> rdm1(const int irdm) const { return rdm1_.at(irdm); }
    std::shared_ptr<const RDM<1>> rdm1_av() const { return rdm1_av_; }

    // returns an occ-occ sized 1RDM
    std::shared_ptr<Matrix> rdm1_mat(std::shared_ptr<const RDM<1>> o) const;
    std::shared_ptr<Matrix> rdm1_mat(const int irdm) const { return rdm1_mat(rdm1_[irdm]); }
    std::shared_ptr<Matrix> rdm1_mat() const { return rdm1_mat(rdm1_av_); }

    std::shared_ptr<const RDM<2>> rdm2(const int irdm) const { return rdm2_.at(irdm); }
    std::shared_ptr<const RDM<2>> rdm2_av() const { return rdm2_av_; }

    std::tuple<std::shared_ptr<RDM<3>>, std::shared_ptr<RDM<4>>> compute_rdm34(const int i) const;

    // function to return a CI vectors from orbital info
    std::shared_ptr<const Dvec> civectors() const;
    std::shared_ptr<Dvec> rdm1deriv(const int istate) const;
    std::shared_ptr<Dvec> rdm2deriv(const int istate) const;
    std::shared_ptr<Dvec> rdm3deriv(const int istate) const;
    std::shared_ptr<Dvec> rdm4deriv(const int istate) const;

    // basis-set projection based on SVD
    virtual std::shared_ptr<Reference> project_coeff(const std::shared_ptr<const Geometry>) const;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Reference)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Reference, T>::value>::type> {
    typedef Reference type;
  };
}

#endif
