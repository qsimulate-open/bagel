//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fock.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;


// Non-DF Fock matrix, standard basis
template <>
void Fock<0>::fock_two_electron_part(shared_ptr<const Matrix> den) {
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

  const int shift = sizeof(int) * 4;
  const int size = basis.size();

  // first make max_density_change vector for each batch pair.
  const double* density_data = density_->data();

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
          cmax = max(cmax, fabs(density_data[iin + jj]));
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
    for (int i1 = i0; i1 != size; ++i1) {
      const unsigned int i01 = i0 *size + i1;

      const shared_ptr<const Shell>  b1 = basis[i1];
      const int b1offset = offset[i1];
      const int b1size = b1->nbasis();

      const double density_change_01 = max_density_change[i01] * 4.0;

      for (int i2 = i0; i2 != size; ++i2) {
        const shared_ptr<const Shell>  b2 = basis[i2];
        const int b2offset = offset[i2];
        const int b2size = b2->nbasis();

        const double density_change_02 = max_density_change[i0 * size + i2];
        const double density_change_12 = max_density_change[i1 * size + i2];

        for (int i3 = i2; i3 != size; ++i3) {
          const unsigned int i23 = i2 * size + i3;
          if (i23 < i01) continue;

          const double density_change_23 = max_density_change[i2 * size + i3] * 4.0;
          const double density_change_03 = max_density_change[i0 * size + i2];
          const double density_change_13 = max_density_change[i0 * size + i2];

          const bool eqli01i23 = (i01 == i23);

          const shared_ptr<const Shell>  b3 = basis[i3];
          const int b3offset = offset[i3];
          const int b3size = b3->nbasis();

          const double mulfactor = max(max(max(density_change_01, density_change_02),
                                           max(density_change_12, density_change_23)),
                                           max(density_change_03, density_change_13));
          const double integral_bound = mulfactor * schwarz_[i01] * schwarz_[i23];
          const bool skip_schwarz = integral_bound < schwarz_thresh_;
          if (skip_schwarz) continue;

          array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
#ifdef LIBINT_INTERFACE
          Libint eribatch(input);
#else
          ERIBatch eribatch(input, mulfactor);
#endif
          eribatch.compute();
          const double* eridata = eribatch.data();
          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int j0n = j0 * ndim();

            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const unsigned int nj01 = (j0 << shift) + j1;
              const bool skipj0j1 = (j0 > j1);
              if (skipj0j1) {
                eridata += b2size * b3size;
                continue;
              }

              const bool eqlj0j1 = (j0 == j1);
              const double scal01 = (eqlj0j1 ? 0.5 : 1.0);
              const int j1n = j1 * ndim();

              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int maxj1j2 = max(j1, j2);
                const int minj1j2 = min(j1, j2);

                const int maxj0j2 = max(j0, j2);
                const int minj0j2 = min(j0, j2);
                const int j2n = j2 * ndim();

                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                  const bool skipj2j3 = (j2 > j3);
                  const unsigned int nj23 = (j2 << shift) + j3;
                  const bool skipj01j23 = (nj01 > nj23) && eqli01i23;

                  if (skipj2j3 || skipj01j23) continue;

                  const int maxj1j3 = max(j1, j3);
                  const int minj1j3 = min(j1, j3);

                  double intval = *eridata * scal01 * (j2 == j3 ? 0.5 : 1.0) * (nj01 == nj23 ? 0.25 : 0.5); // 1/2 in the Hamiltonian absorbed here
                  const double intval4 = 4.0 * intval;

                  element(j1, j0) += density_data[j2n + j3] * intval4;
                  element(j3, j2) += density_data[j0n + j1] * intval4;
                  element(j3, j0) -= density_data[j1n + j2] * intval;
                  element(maxj1j2, minj1j2) -= density_data[j0n + j3] * intval;
                  element(maxj0j2, minj0j2) -= density_data[j1n + j3] * intval;
                  element(maxj1j3, minj1j3) -= density_data[j0n + j2] * intval;
                }
              }
            }
          }

        }
      }
    }
  }
  for (int i = 0; i != ndim(); ++i) element(i, i) *= 2.0;
  fill_upper();
}


template<int DF>
void Fock<DF>::fock_two_electron_part(shared_ptr<const Matrix> den_ex) {
  static_assert(DF == 1, "Wrong Fock matrix implementation is being compiled.");
#ifndef NDEBUG
  cout << "    .. warning .. use a new Fock builder if possible (coeff_ required)" << endl;
#endif

  shared_ptr<const DFDist> df = geom_->df();

  // some constants
  assert(ndim() == df->nbasis0());

  Timer pdebug(3);

  shared_ptr<Matrix> coeff = den_ex->copy();
  *coeff *= -1.0;
  int nocc = 0;
  {
    VectorB vec(ndim());
    coeff->diagonalize(vec);
    for (int i = 0; i != ndim(); ++i) {
      if (vec[i] < -1.0e-8) {
        ++nocc;
        const double fac = std::sqrt(-vec(i));
        for_each(coeff->element_ptr(0,i), coeff->element_ptr(0,i+1), [&fac](double& i) { i *= fac; });
      } else { break; }
    }
  }
  if (nocc == 0) return;
  pdebug.tick_print("Compute coeff (redundant)");

  shared_ptr<DFHalfDist> halfbj = df->compute_half_transform(coeff->slice(0,nocc));
  pdebug.tick_print("First index transform");

  shared_ptr<DFHalfDist> half = halfbj->apply_J();
  pdebug.tick_print("Metric multiply");

  *this += *half->form_2index(half, -0.5);
  pdebug.tick_print("Exchange build");

  *this += *df->compute_Jop(density_);
  pdebug.tick_print("Coulomb build");
}


template<int DF>
void Fock<DF>::fock_two_electron_part_with_coeff(const MatView ocoeff, const bool rhf, const double scale_exchange, const double scale_coulomb) {
  if (DF == 0) throw logic_error("Fock<DF>::fock_two_electron_part_with_coeff() is only for DF cases");

  Timer pdebug(3);

  shared_ptr<const DFDist> df = geom_->df();

  if (scale_exchange != 0.0) {
    shared_ptr<DFHalfDist> halfbj = df->compute_half_transform(ocoeff);
    pdebug.tick_print("First index transform");

    shared_ptr<DFHalfDist> half = halfbj->apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->form_2index(half, -1.0*scale_exchange);
    pdebug.tick_print("Exchange build");

    if (rhf) {
      Matrix oc(ocoeff);
      auto coeff = make_shared<const Matrix>(*oc.transpose()*(2.0*scale_coulomb));
      *this += *df->compute_Jop(half, coeff, true);
    } else {
      shared_ptr<Matrix> jop = df->compute_Jop(density_);
      if (scale_coulomb != 1.0)
        jop->scale(scale_coulomb);
      *this += *jop;
    }
    // when gradient is requested..
    if (store_half_)
      half_ = half;
  } else {
    shared_ptr<Matrix> jop = df->compute_Jop(density_);
    if (scale_coulomb != 1.0)
      jop->scale(scale_coulomb);
    *this += *jop;
  }
  pdebug.tick_print("Coulomb build");

}


template class bagel::Fock<0>;
template class bagel::Fock<1>;

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Fock<0>)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::Fock<1>)

