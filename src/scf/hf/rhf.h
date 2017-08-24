//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rhf.h
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


#ifndef __BAGEL_SRC_SCF_HF_RHF_H
#define __BAGEL_SRC_SCF_HF_RHF_H

#include <src/util/math/diis.h>
#include <src/scf/scf_base.h>
#include <src/scf/levelshift.h>

namespace bagel {

class RHF : public SCF_base {
  protected:
    double lshift_;
    std::shared_ptr<LevelShift<DistMatrix>> levelshift_;

    bool dodf_;
    bool restarted_;

    std::shared_ptr<DIIS<DistMatrix>> diis_;
    std::shared_ptr<const Matrix> compute_Fock_FMM(std::shared_ptr<const Matrix> density, std::shared_ptr<const Matrix> coeff = nullptr);

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<SCF_base>(*this);
      ar << lshift_ << dodf_ << diis_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<SCF_base>(*this);
      ar >> lshift_ >> dodf_ >> diis_;
      if (lshift_ != 0.0)
        levelshift_ = std::make_shared<ShiftVirtual<DistMatrix>>(nocc_, lshift_);
      restarted_ = true;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    RHF() { }
    RHF(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re = nullptr);
    virtual ~RHF() { }

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    bool dodf() const { return dodf_; }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RHF)

#endif
