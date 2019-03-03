//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fock_london.cc
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

#include <src/scf/giaohf/fock_london.h>

using namespace std;
using namespace bagel;


// Non-DF Fock matrix, GIAO basis
template <>
void Fock_London<0>::fock_two_electron_part(shared_ptr<const ZMatrix> den_ex) {
  const vector<shared_ptr<const Atom>> atoms = geom_->atoms();
  vector<shared_ptr<const Shell>> basis;
  vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const vector<shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const vector<int> tmpoff = geom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }
  const int size = basis.size();

  // first make max_density_change vector for each batch pair.
  const complex<double>* density_data = density_->data();

  vector<double> max_density_change(size * size);
  for (int i = 0; i != size; ++i) {
    const int ioffset = offset[i];
    const int isize = basis[i]->nbasis();
    for (int j = i; j != size; ++j) {
      const int joffset = offset[j];
      const int jsize = basis[j]->nbasis();

      double cmax = 0.0;
      for (int ii = ioffset; ii != ioffset + isize; ++ii) {
        const int iin = ii * ndim();
        for (int jj = joffset; jj != joffset + jsize; ++jj) {
          cmax = std::max(cmax, std::abs(density_data[iin + jj]));
        }
      }
      const int ij = i * size + j;
      const int ji = j * size + i;
      max_density_change[ij] = cmax;
      max_density_change[ji] = cmax;
    }
  }

  ////////////////////////////////////////////
  // starting 2-e Fock matrix evaluation!
  ////////////////////////////////////////////
  //////////////// ONLY FOR REFERENCES. //////////////////
  for (int i0 = 0; i0 != size; ++i0) {

    const shared_ptr<const Shell>  b0 = basis[i0];
    const int b0offset = offset[i0];
    const int b0size = b0->nbasis();
    for (int i1 = 0; i1 != size; ++i1) {
      const unsigned int i01 = i0 *size + i1;

      const shared_ptr<const Shell>  b1 = basis[i1];
      const int b1offset = offset[i1];
      const int b1size = b1->nbasis();

      const double density_change_01 = max_density_change[i01] * 4.0;

      for (int i2 = 0; i2 != size; ++i2) {
        const shared_ptr<const Shell>  b2 = basis[i2];
        const int b2offset = offset[i2];
        const int b2size = b2->nbasis();

        const double density_change_02 = max_density_change[i0 * size + i2];
        const double density_change_12 = max_density_change[i1 * size + i2];

        for (int i3 = 0; i3 != size; ++i3) {
          const unsigned int i23 = i2 * size + i3;

          const double density_change_23 = max_density_change[i2 * size + i3] * 4.0;
          const double density_change_03 = max_density_change[i0 * size + i2];
          const double density_change_13 = max_density_change[i0 * size + i2];

          const shared_ptr<const Shell>  b3 = basis[i3];
          const int b3offset = offset[i3];
          const int b3size = b3->nbasis();

          if ((b0offset + b0size + b1offset + b1size) < (b2offset + b3offset)) continue;
          if ((b0offset + b0size + b2offset + b2size) < (b1offset + b3offset)) continue;

          const double mulfactor = std::max(std::max(std::max(density_change_01, density_change_02),
                                           std::max(density_change_12, density_change_23)),
                                           std::max(density_change_03, density_change_13));
          const double integral_bound = mulfactor * schwarz_[i01] * schwarz_[i23];
          const bool skip_schwarz = integral_bound < schwarz_thresh_;
          if (skip_schwarz) continue;

          array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
          ComplexERIBatch eribatch(input, mulfactor);
          eribatch.compute();
          const complex<double>* eridata = eribatch.data();
          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int j0n = j0 * ndim();

            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const int j1n = j1 * ndim();

              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int j2n = j2 * ndim();

                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                  const int j3n = j3 * ndim();
                  complex<double> intval = *eridata * 0.5; // 1/2 in the Hamiltonian absorbed here

                  if (j0 + j1 <  j2 + j3) continue;
                  if (j0 + j2 <  j1 + j3) continue;
                  if (j0 + j1 == j2 + j3) intval *= 0.5;
                  if (j0 + j2 == j1 + j3) intval *= 0.5;
                  const complex<double> intval2 = intval * 2.0;

                  element(j1, j0) += density_data[j3n + j2] * intval2; // Coulomb  (ab|cd)
                  element(j3, j0) -= density_data[j1n + j2] * intval;  // Exchange (ad|cb)

                  element(j3, j2) += density_data[j1n + j0] * intval2; // Coulomb  (cd|ab)
                  element(j1, j2) -= density_data[j3n + j0] * intval;  // Exchange (cb|ad)

                  element(j0, j1) += density_data[j2n + j3] * std::conj(intval2); // Coulomb  (ba|dc)
                  element(j0, j3) -= density_data[j2n + j1] * std::conj(intval);  // Exchange (da|bc)

                  element(j2, j3) += density_data[j0n + j1] * std::conj(intval2); // Coulomb  (dc|ba)
                  element(j2, j1) -= density_data[j0n + j3] * std::conj(intval);  // Exchange (bc|dz)
                }
              }
            }
          }

        }
      }
    }
  }
}


template<int DF>
void Fock_London<DF>::fock_two_electron_part(shared_ptr<const ZMatrix> den_ex) {
  static_assert(DF == 1, "Wrong Fock matrix implementation is being compiled.");
#ifndef NDEBUG
  cout << "    .. warning .. use a new Fock builder if possible (coeff_ required)" << endl;
#endif

  auto df = dynamic_pointer_cast<const ComplexDFDist>(geom_->df());
  assert(df);

  // some constants
  assert(ndim() == df->nbasis0());

  Timer pdebug(3);

  shared_ptr<ZMatrix> coeff = den_ex->copy();
  *coeff *= -1.0;
  int nocc = 0;
  {
    VectorB vec(ndim());
    coeff->diagonalize(vec);
    for (int i = 0; i != ndim(); ++i) {
      if (vec[i] < -1.0e-8) {
        ++nocc;
        const double fac = std::sqrt(-vec(i));
        for_each(coeff->element_ptr(0,i), coeff->element_ptr(0,i+1), [&fac](complex<double>& i) { i *= fac; });
      } else { break; }
    }
  }
  if (nocc == 0) return;
  pdebug.tick_print("Compute coeff (redundant)");

  shared_ptr<ComplexDFHalfDist> halfbj = df->complex_compute_half_transform(coeff->slice(0,nocc));
  pdebug.tick_print("First index transform");

  shared_ptr<ComplexDFHalfDist> half = halfbj->complex_apply_J();
  pdebug.tick_print("Metric multiply");

  *this += *half->complex_form_2index(half, -0.5);
  pdebug.tick_print("Exchange build");

  *this += *df->complex_compute_Jop(density_);
  pdebug.tick_print("Coulomb build");
}


template<int DF>
void Fock_London<DF>::fock_two_electron_part_with_coeff(const ZMatView ocoeff, const bool rhf, const double scale_exchange) {
  if (DF == 0) throw logic_error("Fock_London<DF>::fock_two_electron_part_with_coeff() is only for DF cases");

  Timer pdebug(3);

  auto df = dynamic_pointer_cast<const ComplexDFDist>(geom_->df());
  assert(df);

  if (scale_exchange != 0.0) {
    shared_ptr<ComplexDFHalfDist> halfbj = df->complex_compute_half_transform(ocoeff);
    pdebug.tick_print("First index transform");

    shared_ptr<ComplexDFHalfDist> half = halfbj->complex_apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->complex_form_2index(half, -1.0*scale_exchange);
    pdebug.tick_print("Exchange build");

    if (rhf) {
      auto oc = make_shared<ZMatrix>(ocoeff);
      auto coeff = make_shared<const ZMatrix>(*oc->transpose()*2.0);
      *this += *df->complex_compute_Jop(half, coeff, true);
    } else {
      *this += *df->complex_compute_Jop(density_);
      throw runtime_error("So far, only RHF has been set up with London orbitals, so this should not be called.");
    }
    // when gradient is requested..
    if (store_half_)
      half_ = half;
  } else {
    *this += *df->complex_compute_Jop(density_);
    throw runtime_error("scale_exchange == 0.0 ??  This should not be the case for any London orbital methods.");
  }
  pdebug.tick_print("Coulomb build");
}


template class bagel::Fock_London<0>;
template class bagel::Fock_London<1>;

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Fock_London<0>)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Fock_London<1>)

