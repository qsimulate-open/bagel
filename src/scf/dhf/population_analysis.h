//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: population_analysis.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_SCF_DHF_POPULATION_ANALYSIS_H
#define __SRC_SCF_DHF_POPULATION_ANALYSIS_H

#include <cassert>
#include <src/wfn/geometry.h>
#include <src/util/math/zmatrix.h>
#include <src/util/atommap.h>

namespace bagel {
namespace {

void population_analysis(std::shared_ptr<const Geometry> geom, const ZMatView coeff, std::shared_ptr<const ZMatrix> overlap) {

  const ZMatrix right = *overlap * coeff;
  const ZMatView left = coeff;

  const int nbasis = geom->nbasis();
  const int norb = geom->nbasis();
  assert(nbasis == coeff.ndim()/4);


  for (int i = 0; i < norb; ++i) {

    const int offset = i;
    const int blocksize = 1;

    double alphapart = 0.0;
    double betapart = 0.0;
    double pospart = 0.0;
    double negpart = 0.0;

    const std::complex<double>* left_ptr1 = left.element_ptr(0,i+offset);
    const std::complex<double>* right_ptr1 = right.element_ptr(0,i+offset);
    const std::complex<double>* left_ptr2 = left.element_ptr(0,i+offset+blocksize);
    const std::complex<double>* right_ptr2 = right.element_ptr(0,i+offset+blocksize);

    std::map<std::string, int> element_count;
    int current_ao = 0;
    for (auto& atom : geom->atoms()) {
      std::stringstream base_name;
      base_name << "spatial MO " << i + 1 << " " << atom->name() << "_" << element_count[atom->name()]++;
      for (auto& shell : atom->shells()) {
        if (!shell->dummy()) {
          std::stringstream ss;
          AtomMap am;
          ss << base_name.str() << ":" << am.angular_string(shell->angular_number());

          double val1 = blas::dot_product(left_ptr1+current_ao, shell->nbasis(), right_ptr1+current_ao).real();
          double val2 = blas::dot_product(left_ptr1+current_ao+nbasis, shell->nbasis(), right_ptr1+current_ao+nbasis).real();
          double val3 = blas::dot_product(left_ptr1+current_ao+2*nbasis, shell->nbasis(), right_ptr1+current_ao+2*nbasis).real();
          double val4 = blas::dot_product(left_ptr1+current_ao+3*nbasis, shell->nbasis(), right_ptr1+current_ao+3*nbasis).real();

          // average the alpha and beta spin orbitals
          val1 += blas::dot_product(left_ptr2+current_ao, shell->nbasis(), right_ptr2+current_ao).real();
          val2 += blas::dot_product(left_ptr2+current_ao+nbasis, shell->nbasis(), right_ptr2+current_ao+nbasis).real();
          val3 += blas::dot_product(left_ptr2+current_ao+2*nbasis, shell->nbasis(), right_ptr2+current_ao+2*nbasis).real();
          val4 += blas::dot_product(left_ptr2+current_ao+3*nbasis, shell->nbasis(), right_ptr2+current_ao+3*nbasis).real();
          val1 /= 2.0;
          val2 /= 2.0;
          val3 /= 2.0;
          val4 /= 2.0;

          const double out = val1 + val2 + val3 + val4;
          std::cout << std::setw(8) << ss.str() << std::setw(16) << std::setprecision(8) << out << std::endl;
          current_ao += shell->nbasis();

          alphapart += val1;
          alphapart += val3;
          betapart += val2;
          betapart += val4;
          pospart += val1;
          pospart += val2;
          negpart += val3;
          negpart += val4;
        }
      }
    }
    std::cout << std::endl;
    std::cout << "spatial MO " << i + 1 << " Positive energy basis functions: " << pospart << std::endl;
    std::cout << "spatial MO " << i + 1 << " Negative energy basis functions: " << negpart << std::endl;
    std::cout << std::endl;
  }

}


// Non-relativistic analogue...  Designed for RHF; needs to be adjusted for use w/ UHF, etc.
void population_analysis(std::shared_ptr<const Geometry> geom, const MatView coeff, std::shared_ptr<const Matrix> overlap) {

  const Matrix right = *overlap * coeff;
  const MatView left = coeff;

  const int nbasis = geom->nbasis();
  assert(nbasis == coeff.ndim());

  for (int i = 0; i < nbasis; ++i) {

    double total = 0.0;
    const double* left_ptr1 = left.element_ptr(0,i);
    const double* right_ptr1 = right.element_ptr(0,i);

    std::map<std::string, int> element_count;
    int current_ao = 0;
    for (auto& atom : geom->atoms()) {
      std::stringstream base_name;
      base_name << "spatial MO " << i + 1 << " " << atom->name() << "_" << element_count[atom->name()]++;
      for (auto& shell : atom->shells()) {
        if (!shell->dummy()) {
          std::stringstream ss;
          AtomMap am;
          ss << base_name.str() << ":" << am.angular_string(shell->angular_number());

          double val = blas::dot_product(left_ptr1+current_ao, shell->nbasis(), right_ptr1+current_ao);

          std::cout << std::setw(8) << ss.str() << std::setw(16) << std::setprecision(8) << val << std::endl;
          current_ao += shell->nbasis();
          total += val;
        }
      }
    }
    std::cout << std::endl;
    std::cout << "spatial MO " << i + 1 << " basis functions: " << total << std::endl;
    std::cout << std::endl;
  }

}

}
}

#endif
