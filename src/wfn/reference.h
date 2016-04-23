//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reference.h
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


#ifndef _BAGEL_WFN_REFERENCE_H
#define _BAGEL_WFN_REFERENCE_H

#include <set>
#include <src/scf/coeff.h>
#include <src/mat1e/hcore.h>
#include <src/ci/fci/dvec.h>
#include <src/wfn/geometry.h>
#include <src/wfn/ciwfn.h>
#include <src/wfn/rdm.h>
#include <src/util/vec.h>

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

    std::vector<double> energy_;

    std::shared_ptr<const Hcore> hcore_;
    VectorB eig_;

    int nclosed_;
    int nact_;
    int nvirt_;

    int nstate_;

    std::shared_ptr<const CIWfn> ciwfn_;
    std::shared_ptr<const VecRDM<1>> rdm1_;
    std::shared_ptr<const VecRDM<2>> rdm2_;
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
              const std::vector<double> en,
              std::shared_ptr<const VecRDM<1>> rdm1 = std::make_shared<VecRDM<1>>(),
              std::shared_ptr<const VecRDM<2>> rdm2 = std::make_shared<VecRDM<2>>(),
              std::shared_ptr<const RDM<1>> rdm1_av = nullptr,
              std::shared_ptr<const RDM<2>> rdm2_av = nullptr,
              std::shared_ptr<const CIWfn> ci = nullptr);

    // if only given one energy
    Reference(std::shared_ptr<const Geometry> g, std::shared_ptr<const Coeff> c,
              const int nclo, const int nact, const int nvirt, const double en = 0.0,
              std::shared_ptr<const VecRDM<1>> rdm1 = std::make_shared<VecRDM<1>>(),
              std::shared_ptr<const VecRDM<2>> rdm2 = std::make_shared<VecRDM<2>>(),
              std::shared_ptr<const RDM<1>> rdm1_av = nullptr,
              std::shared_ptr<const RDM<2>> rdm2_av = nullptr,
              std::shared_ptr<const CIWfn> ci = nullptr) :
      Reference(g, c, nclo, nact, nvirt, std::vector<double>(1, en), rdm1, rdm2, rdm1_av, rdm2_av, ci) { }

    // copy construct with optionally updating coeff
    Reference(const Reference& o, std::shared_ptr<const Coeff> c = nullptr) :
      Reference(o.geom(), c ? c : o.coeff(), o.nclosed(), o.nact(), o.nvirt(), o.energy(), o.rdm1(), o.rdm2(), o.rdm1_av(), o.rdm2_av(), o.ciwfn()) { }

    virtual ~Reference() { }

    std::shared_ptr<const Geometry> geom() const { return geom_; }
    const std::vector<double> schwarz() const { return geom_->schwarz(); }
    virtual std::shared_ptr<const Hcore> hcore() const { return hcore_; }
    virtual std::shared_ptr<const Coeff> coeff() const { return coeff_; }

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
    int nstate() const { return energy_.size(); }

    virtual std::shared_ptr<Reference> extract_state(const int istate) const {
      throw std::runtime_error("Reference::extract_state(int istate) has only been implemented for relativistic wavefunctions.");
    }

    // used in UHF
    void set_coeff_AB(const std::shared_ptr<const Coeff> a, const std::shared_ptr<const Coeff> b);
    std::shared_ptr<const Coeff> coeffA() const { return coeffA_; }
    std::shared_ptr<const Coeff> coeffB() const { return coeffB_; }

    double energy() const { return blas::average(energy_); }
    double energy(const int i) const { return energy_[i]; }

    std::shared_ptr<const CIWfn> ciwfn() const { return ciwfn_; }

    std::shared_ptr<const VecRDM<1>> rdm1() const { return rdm1_; }
    std::shared_ptr<const VecRDM<2>> rdm2() const { return rdm2_; }

    std::shared_ptr<const RDM<1>> rdm1(const int ist) const { return rdm1_->at(ist); }
    std::shared_ptr<const RDM<1>> rdm1_av() const { return rdm1_av_; }

    // returns an occ-occ sized 1RDM
    std::shared_ptr<Matrix> rdm1_mat(std::shared_ptr<const RDM<1>> o) const;
    std::shared_ptr<Matrix> rdm1_mat(const int irdm) const { return rdm1_mat(rdm1_->at(irdm)); }
    std::shared_ptr<Matrix> rdm1_mat() const { return rdm1_mat(rdm1_av_); }

    std::shared_ptr<const RDM<2>> rdm2(const int ist) const { return rdm2_->at(ist); }
    std::shared_ptr<const RDM<2>> rdm2_av() const { return rdm2_av_; }

    std::tuple<std::shared_ptr<const RDM<1>>, std::shared_ptr<const RDM<2>>> rdm12(const int ist, const int jst, const bool recompute = false) const;
    std::tuple<std::shared_ptr<const RDM<3>>, std::shared_ptr<const RDM<4>>> rdm34(const int ist, const int jst) const;
    std::tuple<std::shared_ptr<const RDM<3>>, std::shared_ptr<const RDM<3>>>
      rdm34f(const int ist, const int jst, std::shared_ptr<const Matrix> fock) const;

    // function to return a CI vectors from orbital info
    std::shared_ptr<const Dvec> civectors() const;
    std::shared_ptr<Dvec> rdm1deriv(const int istate) const;
    std::shared_ptr<Dvec> rdm2deriv(const int istate) const;
    // 4RDM derivative is precontracted by the Fock matrix
    std::tuple<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
      rdm34deriv(const int istate, std::shared_ptr<const Matrix> fock, const size_t offset, const size_t size) const;

    // basis-set projection based on SVD
    virtual std::shared_ptr<Reference> project_coeff(const std::shared_ptr<const Geometry>, const bool check_geom_change = true) const;

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
