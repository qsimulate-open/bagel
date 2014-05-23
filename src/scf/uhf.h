//
// BAGEL - Parallel electron correlation program.
// Filename: uhf.h
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


#ifndef __BAGEL_SRC_SCF_UHF_H
#define __BAGEL_SRC_SCF_UHF_H

#include <src/scf/scf_base.h>

// I only implement a DF version
//template<int DF>

namespace bagel {

class UHF : public SCF_base {
  protected:
    std::shared_ptr<const Matrix> aodensity_;
    std::shared_ptr<const Matrix> aodensityA_;
    std::shared_ptr<const Matrix> aodensityB_;
    std::shared_ptr<const Coeff> coeffB_;

    std::unique_ptr<double[]> eigB_;
    double* eigB() { return eigB_.get(); }

    void print_S2(const std::string) const;

    void initial_guess();

  public:
    UHF() { }
    UHF(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> re = nullptr)
      : SCF_base(idata, geom, re) {

      std::cout << indent << "*** Open-shell HF ***" << std::endl << std::endl;

      if (re != nullptr) {
        coeff_  = re->coeffA();
        coeffB_ = re->coeffB();
      }
    }
    virtual ~UHF() { }

    std::tuple<std::shared_ptr<const Matrix>,std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>> form_density_uhf() const;

    virtual void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override;
    // return the natural orbital coefficients and nclosed and nact
    std::tuple<std::shared_ptr<Coeff>, int, std::vector<std::shared_ptr<RDM<1>>>> natural_orbitals() const;

    std::shared_ptr<const Matrix> aodensity() const { return aodensity_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::UHF)

#endif
