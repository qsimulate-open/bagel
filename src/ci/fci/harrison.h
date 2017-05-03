//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: harrison.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_FCI_HARRISONZARRABIAN_H
#define __BAGEL_FCI_HARRISONZARRABIAN_H

#include <src/ci/fci/fci.h>
#include <src/ci/fci/space.h>

namespace bagel {

class HarrisonZarrabian : public FCI {

  protected:
    std::shared_ptr<Space_base> space_;

    virtual void const_denom() override;

    virtual std::shared_ptr<Dvec> form_sigma(std::shared_ptr<const Dvec> c, std::shared_ptr<const MOFile> jop, const std::vector<int>& conv) const override;

    // run-time functions.
    void sigma_aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_bb(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, std::shared_ptr<const MOFile> jop) const;
    void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d) const;
    void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, std::shared_ptr<const MOFile> jop) const;
    void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e) const;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<FCI>(*this);
      std::shared_ptr<const Matrix> coeff = jop_->coeff();
      ar << space_ << coeff;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<FCI>(*this);
      std::shared_ptr<const Matrix> coeff;
      ar >> space_ >> coeff;
      update(coeff);
    }

  public:
    HarrisonZarrabian() { }

    HarrisonZarrabian(std::shared_ptr<const PTree> a, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> b,
        const int ncore = -1, const int nocc = -1, const int nstate = -1, const bool store = false);

    virtual void update(std::shared_ptr<const Matrix>) override;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::HarrisonZarrabian)

#endif

