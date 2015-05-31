//
// BAGEL - Parallel electron correlation program.
// Filename: relreference.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_WFN_RELREFERENCE_H
#define __SRC_WFN_RELREFERENCE_H

#include <src/wfn/reference.h>
#include <src/wfn/relcoeff.h>
#include <src/util/kramers.h>

namespace bagel {

class RelReference : public Reference {
  protected:
    bool gaunt_;
    bool breit_;
    int nneg_;
    std::shared_ptr<const RelCoeff_Striped> relcoeff_;
    bool kramers_;  // Indicates whether or not relcoeff_ has been kramers-adapted

    // RDM things
    std::shared_ptr<const ZMatrix> rdm1_av_;
    std::shared_ptr<const ZMatrix> rdm2_av_;

    std::shared_ptr<const RelCIWfn> ciwfn_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Reference>(*this) & gaunt_ & breit_ & nneg_ & relcoeff_ & kramers_;
    }

  public:
    RelReference() { }
    RelReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const RelCoeff_Striped> c, const double en,
                 const int nneg, const int nocc, const int nact, const int nvirt, const bool ga, const bool br, const bool kram = false,
//               std::shared_ptr<const VecRDM<1>> rdm1 = std::make_shared<VecRDM<1>>(),
//               std::shared_ptr<const VecRDM<2>> rdm2 = std::make_shared<VecRDM<2>>(),
                 std::shared_ptr<const ZMatrix> rdm1_av = nullptr,
                 std::shared_ptr<const ZMatrix> rdm2_av = nullptr,
                 std::shared_ptr<const RelCIWfn> ci = nullptr)
     : Reference(g, nullptr, nocc, nact, nvirt, en), gaunt_(ga), breit_(br), nneg_(nneg), relcoeff_(c), kramers_(kram),
                                                     rdm1_av_(rdm1_av), rdm2_av_(rdm2_av), ciwfn_(ci) {
    }

    std::shared_ptr<const Coeff> coeff() const override { throw std::logic_error("RelReference::coeff() should not be called"); }
    std::shared_ptr<const RelCoeff_Striped> relcoeff() const { return relcoeff_->electronic_part(); }
    std::shared_ptr<const RelCoeff_Striped> relcoeff_full() const { return relcoeff_; }

    bool gaunt() const { return gaunt_; }
    bool breit() const { return breit_; }
    int nneg() const { return nneg_; }
    bool kramers() const { return kramers_; }

    std::shared_ptr<const RelCIWfn> ciwfn() const { return ciwfn_; }

    std::shared_ptr<const ZMatrix> rdm1_av() const { return rdm1_av_; }
    std::shared_ptr<const ZMatrix> rdm2_av() const { return rdm2_av_; }

    std::shared_ptr<const Kramers<2,ZRDM<1>>> rdm1(const int ist, const int jst) const;
    std::shared_ptr<const Kramers<4,ZRDM<2>>> rdm2(const int ist, const int jst) const;
    std::shared_ptr<const Kramers<6,ZRDM<3>>> rdm3(const int ist, const int jst) const;
    std::shared_ptr<const Kramers<8,ZRDM<4>>> rdm4(const int ist, const int jst) const;
    std::shared_ptr<const Kramers<6,ZRDM<3>>> frdm4(const int ist, const int jst, std::shared_ptr<const ZMatrix> fock) const;

    std::shared_ptr<Reference> project_coeff(std::shared_ptr<const Geometry> geomin, const bool check_geom_change = true) const override;

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelReference)

#endif
