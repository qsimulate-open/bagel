//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: shell_ecp.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_MOLECULE_SHELL_ECP_H
#define __SRC_MOLECULE_SHELL_ECP_H

#include <src/molecule/shell_base.h>

namespace bagel {

class Shell_ECP : public Shell_base {

  protected:
    std::vector<double> ecp_exponents_;
    std::vector<double> ecp_coefficients_;
    std::vector<int> ecp_r_power_;

  private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Shell_base>(*this) & ecp_exponents_ & ecp_coefficients_ & ecp_r_power_;
    }

  public:
    Shell_ECP();

    Shell_ECP(const std::array<double,3>& position, const int angular_num, const std::vector<double>& ecp_exponents,
              const std::vector<double>& ecp_coefficients, const std::vector<int>& ecp_r_power);

    double ecp_exponents(const int i) const { return ecp_exponents_[i]; }
    const std::vector<double>& ecp_exponents() const { return ecp_exponents_; }
    const double* ecp_exponents_pointer() const { return &(ecp_exponents_[0]); }

    double ecp_coefficients(const int i) const { return ecp_coefficients_[i]; }
    const std::vector<double>& ecp_coefficients() const { return ecp_coefficients_; }
    const double* ecp_coefficients_pointer() const { return &(ecp_coefficients_[0]); }

    int ecp_r_power(const int i) const { return ecp_r_power_[i]; }
    const std::vector<int>& ecp_r_power() const { return ecp_r_power_; }
    const int* ecp_r_power_pointer() const { return &(ecp_r_power_[0]); }

    std::string show() const override;

};

}

#endif

