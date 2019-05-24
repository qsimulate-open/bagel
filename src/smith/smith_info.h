//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith_info.h
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

#ifndef __SRC_SMITH_SMITH_INFO_H
#define __SRC_SMITH_SMITH_INFO_H

#include <src/wfn/relreference.h>

namespace bagel {

template<typename DataType>
class SMITH_Info {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
    template<int N>
    using RDMType = typename std::conditional<std::is_same<DataType,double>::value,RDM<N>,Kramers<N*2,ZRDM<N>>>::type;
    using CIWfnT  = typename std::conditional<std::is_same<DataType,double>::value,CIWfn,RelCIWfn>::type;

    std::shared_ptr<const Reference> ref_;
    std::string method_;

    int ncore_;
    int nfrozenvirt_;
    double thresh_;
    double shift_;
    int maxiter_;
    int maxtile_;
    size_t cimaxchunk_;
    int davidson_subspace_;

    bool grad_;

    bool do_ms_;
    bool do_xms_;
    bool sssr_;
    bool shift_diag_;
    bool shift_imag_;
    bool block_diag_fock_;
    // use orthogonal basis in linear equation for CASPT2 case
    bool orthogonal_basis_;

    bool restart_;
    bool restart_each_iter_;
    bool convergence_throw_;

    double thresh_overlap_;

    // For restarted jobs
    int state_begin_;
    int restart_iter_;

    std::shared_ptr<const PTree> aniso_data_;  // Inputs to pseudospin Hamiltonian module
    std::string external_rdm_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & ref_ & method_ & ncore_ & nfrozenvirt_ & thresh_ & shift_ & maxiter_;
      ar & maxtile_ & cimaxchunk_ & davidson_subspace_ & grad_;
      ar & do_ms_ & do_xms_ & sssr_ & shift_diag_ & shift_imag_ & block_diag_fock_ & orthogonal_basis_ & restart_ & restart_each_iter_ & convergence_throw_;
      ar & thresh_overlap_ & state_begin_ & restart_iter_ & aniso_data_ & external_rdm_;
    }

  public:
    SMITH_Info() { }
    SMITH_Info(std::shared_ptr<const Reference> o, const std::shared_ptr<const PTree> idata);
    SMITH_Info(std::shared_ptr<const Reference> o, std::shared_ptr<const SMITH_Info> info);

    std::string method() const { return method_; }
    int ncore() const { return ncore_; }
    int nclosed() const { return ref_->nclosed(); }
    int nact() const { return ref_->nact(); }
    int nocc() const { return ref_->nocc(); }
    int nvirt() const { return ref_->nvirt() - nfrozenvirt_; }
    int nfrozenvirt() const { return nfrozenvirt_; }

    std::shared_ptr<const MatType> hcore() const { return nullptr; }

    std::shared_ptr<const RDM<1,DataType>> rdm1_av() const;

    std::tuple<std::shared_ptr<const RDMType<1>>, std::shared_ptr<const RDMType<2>>> rdm12(const int ist, const int jst) const;
    std::tuple<std::shared_ptr<const RDMType<3>>, std::shared_ptr<RDMType<3>>> rdm34f(const int ist, const int jst, std::shared_ptr<const MatType> fock) const;
    std::tuple<std::shared_ptr<const RDMType<3>>, std::shared_ptr<const RDMType<4>>> rdm34(const int ist, const int jst) const;
    std::shared_ptr<RDMType<3>> rdm4f_contract(std::shared_ptr<const RDMType<3>>, std::shared_ptr<const RDMType<4>>, std::shared_ptr<const MatType> fock) const;

    double thresh() const { return thresh_; }
    double shift() const {return shift_; }
    int maxiter() const { return maxiter_; }
    int maxtile() const { return maxtile_; }
    int cimaxchunk() const { return cimaxchunk_; }
    bool grad() const { return grad_; }
    bool do_ms() const { return do_ms_; }
    bool do_xms() const { return do_xms_; }
    bool sssr() const { return sssr_; }
    bool shift_diag() const { return shift_diag_; }
    bool shift_imag() const { return shift_imag_; }
    bool block_diag_fock() const { return block_diag_fock_; }
    bool restart() const { return restart_; }
    bool restart_each_iter() const { return restart_each_iter_; }
    bool convergence_throw() const { return convergence_throw_; }
    bool orthogonal_basis() const { return orthogonal_basis_; }
    bool rdm4_eval() const { return (grad_ || method_!="caspt2"); }

    double thresh_overlap() const { return thresh_overlap_; }

    template<typename T = DataType, class = typename std::enable_if<std::is_same<T, std::complex<double>>::value>::type>
    bool gaunt() const { return relref()->gaunt(); }
    template<typename T = DataType, class = typename std::enable_if<std::is_same<T, std::complex<double>>::value>::type>
    bool breit() const { return relref()->breit(); }

    template<typename T = DataType, class = typename std::enable_if<std::is_same<T, std::complex<double>>::value>::type>
    std::shared_ptr<const RelReference> relref() const {
      assert(std::dynamic_pointer_cast<const RelReference>(ref_));
      return std::dynamic_pointer_cast<const RelReference>(ref_);
    }

    int davidson_subspace() const { return davidson_subspace_; }

    std::shared_ptr<const Reference> ref() const { return ref_; }
    std::shared_ptr<const Geometry> geom() const { return ref_->geom(); }
    std::shared_ptr<const CIWfnT> ciwfn() const;
    std::string external_rdm() const { return external_rdm_; }

    // this function hides coeff function in Reference and RelReference
    std::shared_ptr<const MatType> coeff() const { assert(false); }

    // TODO Do we want to keep this?  Implemented for debugging, but could be useful in the future
    std::shared_ptr<const Reference> extract_ref(const std::vector<int> states, const bool extract_rdm) const;

    int state_begin() const { return state_begin_; }
    int restart_iter() const { return restart_iter_; }

    void set_restart_params(const int state, const int iter) {
      state_begin_ = state;
      restart_iter_ = iter;
      if (state_begin_ < 0 || state_begin_ > (nact() ? ciwfn()->nstates() : 1) || restart_iter_ < 0)
        throw std::runtime_error("Invalid starting point for RelSMITH continue");
    }

    std::shared_ptr<const PTree> aniso_data() const { return aniso_data_; }
};

template<> std::tuple<std::shared_ptr<const RDM<1>>, std::shared_ptr<const RDM<2>>> SMITH_Info<double>::rdm12(const int ist, const int jst) const;
template<> std::tuple<std::shared_ptr<const RDM<3>>, std::shared_ptr<RDM<3>>> SMITH_Info<double>::rdm34f(const int ist, const int jst, std::shared_ptr<const Matrix> fock) const;
template<> std::tuple<std::shared_ptr<const RDM<3>>, std::shared_ptr<const RDM<4>>> SMITH_Info<double>::rdm34(const int ist, const int jst) const;
template<> std::shared_ptr<RDM<3>> SMITH_Info<double>::rdm4f_contract(std::shared_ptr<const RDM<3>> rdm3, std::shared_ptr<const RDM<4>> rdm4, std::shared_ptr<const Matrix> fock) const;
template<> std::tuple<std::shared_ptr<const Kramers<2,ZRDM<1>>>, std::shared_ptr<const Kramers<4,ZRDM<2>>>>
           SMITH_Info<std::complex<double>>::rdm12(const int ist, const int jst) const;
template<> std::tuple<std::shared_ptr<const Kramers<6,ZRDM<3>>>, std::shared_ptr<Kramers<6,ZRDM<3>>>>
           SMITH_Info<std::complex<double>>::rdm34f(const int ist, const int jst, std::shared_ptr<const ZMatrix> fock) const;
template<> std::tuple<std::shared_ptr<const Kramers<6,ZRDM<3>>>, std::shared_ptr<const Kramers<8,ZRDM<4>>>>
           SMITH_Info<std::complex<double>>::rdm34(const int ist, const int jst) const;
template<> std::shared_ptr<Kramers<6,ZRDM<3>>> SMITH_Info<std::complex<double>>::rdm4f_contract(std::shared_ptr<const Kramers<6,ZRDM<3>>> rdm3, std::shared_ptr<const Kramers<8,ZRDM<4>>> rdm4, std::shared_ptr<const ZMatrix> fock) const;

template<> std::shared_ptr<const CIWfn>   SMITH_Info<double>::ciwfn() const;
template<> std::shared_ptr<const Matrix>  SMITH_Info<double>::coeff() const;
template<> std::shared_ptr<const Matrix>  SMITH_Info<double>::hcore() const;
template<> std::shared_ptr<const RDM<1>>  SMITH_Info<double>::rdm1_av() const;
template<> std::shared_ptr<const RelCIWfn>SMITH_Info<std::complex<double>>::ciwfn() const;
template<> std::shared_ptr<const ZMatrix> SMITH_Info<std::complex<double>>::coeff() const;
template<> std::shared_ptr<const ZMatrix> SMITH_Info<std::complex<double>>::hcore() const;
template<> std::shared_ptr<const ZRDM<1>> SMITH_Info<std::complex<double>>::rdm1_av() const;

extern template class SMITH_Info<double>;
extern template class SMITH_Info<std::complex<double>>;

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SMITH_Info<double>)
BOOST_CLASS_EXPORT_KEY(bagel::SMITH_Info<std::complex<double>>)

#endif
