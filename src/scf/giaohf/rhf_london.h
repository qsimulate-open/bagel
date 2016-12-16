//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rhf_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_LONDON_RHF_LONDON_H
#define __BAGEL_SRC_LONDON_RHF_LONDON_H

#include <src/scf/scf_base.h>
#include <src/scf/levelshift.h>
#include <src/scf/giaohf/fock_london.h>
#include <src/util/math/diis.h>
#include <src/wfn/method.h>

namespace bagel {

class RHF_London : public SCF_base_London {
  protected:
    double lshift_;
    std::shared_ptr<LevelShift<DistZMatrix>> levelshift_;

    bool dodf_;
    bool restarted_;

    std::shared_ptr<DIIS<DistZMatrix,ZMatrix>> diis_;

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<SCF_base_London>(*this);
      ar << lshift_ << dodf_ << diis_;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<SCF_base_London>(*this);
      ar >> lshift_ >> dodf_ >> diis_;
      if (lshift_ != 0.0)
        levelshift_ = std::make_shared<ShiftVirtual<DistZMatrix>>(nocc_, lshift_);
      restarted_ = true;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

  public:
    RHF_London() { }
    RHF_London(const std::shared_ptr<const PTree> idata_, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re = nullptr);
    virtual ~RHF_London() { }

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;

    bool dodf() const { return dodf_; }

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RHF_London)

#endif
