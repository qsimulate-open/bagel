//
// BAGEL - Parallel electron correlation program.
// Filename: smith_info.h
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
    int maxiter_;
    int target_;
    int maxtile_;
    int davidson_subspace_;

    bool grad_;

  public:
    SMITH_Info(std::shared_ptr<const Reference> o, const std::shared_ptr<const PTree> idata) : ref_(o) {
      method_ = idata->get<std::string>("method");

      const bool frozen = idata->get<bool>("frozen", true);
      ncore_ = idata->get<int>("ncore", (frozen ? ref_->geom()->num_count_ncore_only()/2 : 0));
      if (ncore_)
        std::cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << std::endl;
      nfrozenvirt_ = idata->get<int>("nfrozenvirt", 0);
      if (nfrozenvirt_)
        std::cout << "    * freezing " << nfrozenvirt_ << " orbital" << (nfrozenvirt_^1 ? "s" : "") << " (virtual)" << std::endl;

      maxiter_ = idata->get<int>("maxiter", 50);
      target_  = idata->get<int>("target",   0);
      maxtile_ = idata->get<int>("maxtile", 10);
      grad_    = idata->get<bool>("grad", false);

      thresh_ = idata->get<double>("thresh", grad_ ? 1.0e-8 : 1.0e-6);
      davidson_subspace_ = idata->get<int>("davidson_subspace", 10);
    }

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
    std::tuple<std::shared_ptr<const RDMType<3>>, std::shared_ptr<const RDMType<4>>> rdm34(const int ist, const int jst) const;
    std::shared_ptr<const RDMType<3>> frdm4(const int ist, const int jst, std::shared_ptr<const MatType> fock) const;

    double thresh() const { return thresh_; }
    int maxiter() const { return maxiter_; }
    int target() const { return target_; }
    int maxtile() const { return maxtile_; }
    bool grad() const { return grad_; }

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

    // this function hides coeff function in Reference and RelReference
    std::shared_ptr<const MatType> coeff() const { assert(false); }
};

template<> std::tuple<std::shared_ptr<const RDM<1>>, std::shared_ptr<const RDM<2>>> SMITH_Info<double>::rdm12(const int ist, const int jst) const;
template<> std::tuple<std::shared_ptr<const RDM<3>>, std::shared_ptr<const RDM<4>>> SMITH_Info<double>::rdm34(const int ist, const int jst) const;
template<> std::tuple<std::shared_ptr<const Kramers<2,ZRDM<1>>>, std::shared_ptr<const Kramers<4,ZRDM<2>>>>
           SMITH_Info<std::complex<double>>::rdm12(const int ist, const int jst) const;
template<> std::tuple<std::shared_ptr<const Kramers<6,ZRDM<3>>>, std::shared_ptr<const Kramers<8,ZRDM<4>>>>
           SMITH_Info<std::complex<double>>::rdm34(const int ist, const int jst) const;

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

#endif
