//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moldenin.h
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

#ifndef __SRC_IO_MOLDENIN_H
#define __SRC_IO_MOLDENIN_H

#include <src/util/io/moldenio.h>
#include <src/wfn/zcoeff.h>

namespace bagel {

class MOInfo {
  public:
    std::shared_ptr<const Geometry> geom;
    MOInfo(std::shared_ptr<const Geometry> g) : geom(g) { }
    std::shared_ptr<Coeff> coeff;
    std::shared_ptr<Coeff> coeffB;
    std::shared_ptr<ZCoeff> zcoeff;
    std::shared_ptr<ZCoeff_Striped> relcoeff;
    VectorB eig;
    VectorB eigB;
    VectorB occup;
    VectorB occupB;
    int norb() const { return coeff ? coeff->mdim() : (zcoeff ? zcoeff->mdim() : relcoeff->mdim()/2); }
    bool has_beta() const { return !eigB.empty(); }
};


class MoldenIn : public MoldenIO {
  protected:
    bool is_spherical_;
    bool cartesian_;

    std::vector<std::shared_ptr<const Atom>> atoms_;

    std::vector<double> mo_eig_;
    std::vector<double> mo_occup_;
    std::vector<int> mo_spin_;
    std::vector<std::vector<double>> mo_coefficients_;
    std::vector<std::vector<std::complex<double>>> mo_coefficients_complex_;
    std::vector<std::vector<molden_impl::complex4>> mo_coefficients_relativistic_;

    std::vector<std::vector<std::vector<std::pair<int, double>>>> lmtuv_;
    std::vector<int> gto_order_;
    std::vector<std::vector<int>> shell_orders_;

    void compute_transforms();
    void read_mos(MOInfo&);
    void read_mos_complex(MOInfo&);
    void read_mos_relativistic(MOInfo&);

  public:
    MoldenIn(const std::string filename, const bool is_spherical = true);

    void read();

    MoldenIn& operator>> (std::vector<std::shared_ptr<const Atom>>& atoms_);
    MoldenIn& operator>> (MOInfo& mo);
    bool has_mo() const { return !mo_coefficients_.empty() || !mo_coefficients_complex_.empty() || !mo_coefficients_relativistic_.empty(); }
    bool has_beta() const { return std::any_of(mo_spin_.begin(), mo_spin_.end(), [](const int i){return i==1;}); }

  protected:
    template<typename Type>
    void moread_util(const Type* mo, Type* idata, const std::vector<int>& atom_offsets) {
      int ii = 0;
      for (auto iatom = shell_orders_.begin(); iatom != shell_orders_.end(); ++iatom, ++ii) {
        if (iatom->empty()) continue;
        Type* tmp_idata = idata + atom_offsets[gto_order_[ii]-1];
        for (auto& ishell : *iatom) {
          if (cartesian_) {
            std::vector<int> corder = m2b_cart_.at(ishell);
            std::vector<double> scales = scaling_.at(ishell);
            int jj = 0;
            if (is_spherical_) {
              std::vector<Type> in;
              for(auto& iorder : corder)
                in.push_back(*(mo + iorder) * scales.at(jj++));
              std::vector<Type> new_in = transform_cart(in, ishell);
              tmp_idata = std::copy(new_in.begin(), new_in.end(), tmp_idata);
            } else {
              for (auto& iorder : corder)
                *tmp_idata++ = *(mo + iorder) * scales.at(jj++);
            }
            mo += corder.size();
          } else {
            std::vector<int> corder = m2b_sph_.at(ishell);
            for (auto& iorder : corder)
              *tmp_idata++ = *(mo + iorder);
            mo += corder.size();
          }
        }
      }
    }

    template<typename Type>
    std::vector<Type> transform_cart(std::vector<Type> carts, int ang_l) const {
      std::vector<std::vector<std::pair<int,double>>> mtuv = lmtuv_.at(ang_l);
      std::vector<Type> out;
      for (auto& im : mtuv) {
        Type value = 0.0;
        for (auto& ituv : im)
          value += carts.at(ituv.first) * ituv.second;
        out.push_back(value);
      }
      return out;
    }
};

}

#endif
