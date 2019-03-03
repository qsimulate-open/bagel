//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenout.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __io_moldenout_h
#define __io_moldenout_h

#include <src/util/io/moldenio.h>

namespace bagel {

class MoldenOut : public MoldenIO {
   protected:
      std::ofstream ofs_;

      void write_geom();
      void write_aos();
      void write_mos();
      void write_mos_complex();
      void write_mos_relativistic();
      void write_freq();

      template<typename Type>
      void write_mo_one(std::ostream& ss, const Type* data) {
        const bool is_spherical = mol_->spherical();
        const int width = std::is_same<Type,double>::value ? 22 : 44;
        int j = 0;
        for (auto& iatom : mol_->atoms()) {
          for (auto& ishell : iatom->shells()) {
            for (int icont = 0; icont != ishell->num_contracted(); ++icont) {
              std::vector<int> corder = (is_spherical ? b2m_sph_.at(ishell->angular_number()) : b2m_cart_.at(ishell->angular_number()));
              std::vector<double> scales = (is_spherical ? std::vector<double>(corder.size(), 1.0) : scaling_.at(ishell->angular_number()));
              for (auto& iorder : corder) {
                ss << std::fixed << std::setw(4) << ++j << " " << std::setw(width) << std::setprecision(16) << data[iorder] / scales.at(iorder) << std::endl;
              }
              data += corder.size();
            }
          }
        }
      }

   public:
      MoldenOut(std::string filename);

      MoldenOut& operator<< (std::shared_ptr<const Molecule>);
      MoldenOut& operator<< (std::shared_ptr<const Reference>);

};

}

#endif
