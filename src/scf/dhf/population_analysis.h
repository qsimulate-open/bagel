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

// Helper function to perform population analysis on 4-component molecular orbitals
// When "paired" is true, we work with spatial orbitals, assuming comes in striped format (1a, 1b, 2a, 2b...)
// When working with spatial MOs, we average the contributions of the two spin-MOs, although they are normally identical

#include <cassert>
#include <src/wfn/geometry.h>
#include <src/util/math/zmatrix.h>
#include <src/util/atommap.h>

namespace bagel {
namespace {

void population_analysis(std::shared_ptr<const Geometry> geom, const ZMatView coeff, std::shared_ptr<const ZMatrix> overlap,
                         const bool paired, const int nclosed = 0, const int nact = 0) {

  const ZMatrix right = *overlap * coeff;
  const ZMatView left = coeff;

  const int nbasis = geom->nbasis();
  const int norb = (!paired) ? 2*geom->nbasis() : geom->nbasis();
  assert(nbasis == coeff.ndim()/4);


  for (int i = 0; i < norb; ++i) {

    const int offset = paired ? i : 0;
    const int blocksize = paired ? 1 : 0;

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
      base_name << (paired ? "spatial " : "spin-") << "MO " << i + 1 << " " << atom->name() << "_" << element_count[atom->name()]++;
      for (auto& shell : atom->shells()) {
        if (!shell->dummy()) {
          std::stringstream ss;
          AtomMap am;
          ss << base_name.str() << ":" << am.angular_string(shell->angular_number());

          double val1 = blas::dot_product(left_ptr1+current_ao, shell->nbasis(), right_ptr1+current_ao).real();
          double val2 = blas::dot_product(left_ptr1+current_ao+nbasis, shell->nbasis(), right_ptr1+current_ao+nbasis).real();
          double val3 = blas::dot_product(left_ptr1+current_ao+2*nbasis, shell->nbasis(), right_ptr1+current_ao+2*nbasis).real();
          double val4 = blas::dot_product(left_ptr1+current_ao+3*nbasis, shell->nbasis(), right_ptr1+current_ao+3*nbasis).real();

          if (paired) {
            // average the alpha and beta spin orbitals
            val1 += blas::dot_product(left_ptr2+current_ao, shell->nbasis(), right_ptr2+current_ao).real();
            val2 += blas::dot_product(left_ptr2+current_ao+nbasis, shell->nbasis(), right_ptr2+current_ao+nbasis).real();
            val3 += blas::dot_product(left_ptr2+current_ao+2*nbasis, shell->nbasis(), right_ptr2+current_ao+2*nbasis).real();
            val4 += blas::dot_product(left_ptr2+current_ao+3*nbasis, shell->nbasis(), right_ptr2+current_ao+3*nbasis).real();
            val1 /= 2.0;
            val2 /= 2.0;
            val3 /= 2.0;
            val4 /= 2.0;
          }

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
    if (!paired) {
      std::cout << (paired ? "spatial " : "spin-") << "MO " << i + 1 << " Alpha basis functions: " << alphapart << std::endl;
      std::cout << (paired ? "spatial " : "spin-") << "MO " << i + 1 << " Beta  basis functions: " << betapart << std::endl;
    }
    std::cout << (paired ? "spatial " : "spin-") << "MO " << i + 1 << " Positive energy basis functions: " << pospart << std::endl;
    std::cout << (paired ? "spatial " : "spin-") << "MO " << i + 1 << " Negative energy basis functions: " << negpart << std::endl;
    std::cout << std::endl;
  }

}

}
}

#endif
