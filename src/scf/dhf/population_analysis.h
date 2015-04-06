//
// BAGEL - Parallel electron correlation program.
// Filename: population_analysis.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_SCF_DHF_POPULATION_ANALYSIS_H
#define __SRC_SCF_DHF_POPULATION_ANALYSIS_H

#include<src/wfn/geometry.h>
#include<src/util/math/zmatrix.h>
#include<src/util/atommap.h>

namespace bagel {
namespace {

void population_analysis(std::shared_ptr<const Geometry> geom, const ZMatView coeff_pos, std::shared_ptr<const ZMatrix> overlap) {

  ZMatrix right = *overlap * coeff_pos;
  const ZMatView left = coeff_pos;

  const int nbasis = geom->nbasis();
  const int norb = coeff_pos.mdim()/2;
  assert(nbasis == coeff_pos.ndim()/4);
  for (int i = 0; i < norb; ++i) {
    const std::complex<double>* left_ptr = left.element_ptr(0,i*2);
    const std::complex<double>* right_ptr = right.element_ptr(0,i*2);

    std::map<std::string, int> element_count;
    int current_ao = 0;
    for (auto& atom : geom->atoms()) {
      std::stringstream base_name;
      base_name << "MO " << i + 1 << " " << atom->name() << "_" << element_count[atom->name()]++;
      for (auto& shell : atom->shells()) {
        if (!shell->dummy()) {
          std::stringstream ss;
          AtomMap am;
          ss << base_name.str() << ":" << am.angular_string(shell->angular_number());
          const double val = blas::dot_product(left_ptr+current_ao, shell->nbasis(), right_ptr+current_ao).real()
                            + blas::dot_product(left_ptr+current_ao+nbasis, shell->nbasis(), right_ptr+current_ao+nbasis).real();

                            //  To include the negative energy contributions as well
                            //+ blas::dot_product(left_ptr+current_ao+2*nbasis, shell->nbasis(), right_ptr+current_ao+2*nbasis).real()
                            //+ blas::dot_product(left_ptr+current_ao+3*nbasis, shell->nbasis(), right_ptr+current_ao+3*nbasis).real();
          current_ao += shell->nbasis();
          std::cout << std::setw(8) << ss.str() << std::setw(16) << std::setprecision(8) << val << std::endl;
        }
      }
    }
    std::cout << std::endl;
  }
}

}
}

#endif
