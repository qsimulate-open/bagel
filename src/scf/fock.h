//
// BAGEL - Parallel electron correlation program.
// Filename: fock.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __BAGEL_SRC_SCF_FOCK_H
#define __BAGEL_SRC_SCF_FOCK_H

#include <src/df/df.h>
#include <src/integral/libint/libint.h>
#include <src/integral/rys/eribatch.h>
#include <src/scf/fock_base.h>

namespace bagel {

template<int DF>
class Fock : public Fock_base {
  protected:
    void fock_two_electron_part(std::shared_ptr<const Matrix> den = std::shared_ptr<Matrix>());
    void fock_two_electron_part_with_coeff(const std::shared_ptr<const Matrix> coeff, const bool rhf, const double scale_ex);

  public:
    // Fock operator for DF cases
    Fock(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const Matrix> b, const std::shared_ptr<const Matrix> c,
         const std::shared_ptr<const Matrix> ocoeff, const bool rhf = false, const double scale_ex = 1.0)
     : Fock_base(a,b,c) {
      fock_two_electron_part_with_coeff(ocoeff, rhf, scale_ex);
      fock_one_electron_part();
    }

    // Fock operator
    Fock(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const Matrix> b, const std::shared_ptr<const Matrix> c, const std::vector<double>& d) : Fock(a,b,c,c,d) {}

    // Fock operator with a different density matrix for exchange
    Fock(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const Matrix> b, const std::shared_ptr<const Matrix> c, std::shared_ptr<const Matrix> ex,
         const std::vector<double>& d)
     : Fock_base(a,b,c,d) {
      fock_two_electron_part(ex);
      fock_one_electron_part();
    }

//  Fock(const std::shared_ptr<const Geometry> a) : Fock_base(a) {}

};


template<int DF>
void Fock<DF>::fock_two_electron_part(std::shared_ptr<const Matrix> den_ex) {

  if (den_ex != density_) throw std::logic_error("den_ex in Fock<DF>::fock_two_electron_part is only with DF");
  if (den_ex == nullptr) den_ex = density_;

  const std::vector<std::shared_ptr<const Atom>> atoms = geom_->atoms();
  std::vector<std::shared_ptr<const Shell>> basis;
  std::vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const std::vector<std::shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const std::vector<int> tmpoff = geom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }

  const int shift = sizeof(int) * 4;
  const int size = basis.size();

  // first make max_density_change std::vector for each batch pair.
  const double* density_data = density_->data();

  std::vector<double> max_density_change(size * size);
  for (int i = 0; i != size; ++i) {
    const int ioffset = offset[i];
    const int isize = basis[i]->nbasis();
    for (int j = i; j != size; ++j) {
      const int joffset = offset[j];
      const int jsize = basis[j]->nbasis();

      double cmax = 0.0;
      for (int ii = ioffset; ii != ioffset + isize; ++ii) {
        const int iin = ii * ndim_;
        for (int jj = joffset; jj != joffset + jsize; ++jj) {
          cmax = std::max(cmax, ::fabs(density_data[iin + jj]));
        }
      }
      const int ij = i * size + j;
      const int ji = i * size + j;
      max_density_change[ij] = cmax;
      max_density_change[ji] = cmax;
    }
  }

  ////////////////////////////////////////////
  // starting 2-e Fock matrix evaluation!
  ////////////////////////////////////////////
  if (DF == 0) {
    //////////////// ONLY FOR REFERENCES. //////////////////
    std::shared_ptr<Petite> plist = geom_->plist();;
    const bool c1 = plist->nirrep() == 1;

    for (int i0 = 0; i0 != size; ++i0) {
      if (!plist->in_p1(i0)) continue;

      const std::shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = offset[i0];
      const int b0size = b0->nbasis();
      for (int i1 = i0; i1 != size; ++i1) {
        const unsigned int i01 = i0 *size + i1;
        if (!plist->in_p2(i01)) continue;

        const std::shared_ptr<const Shell>  b1 = basis[i1];
        const int b1offset = offset[i1];
        const int b1size = b1->nbasis();

        const double density_change_01 = max_density_change[i01] * 4.0;

        for (int i2 = i0; i2 != size; ++i2) {
          const std::shared_ptr<const Shell>  b2 = basis[i2];
          const int b2offset = offset[i2];
          const int b2size = b2->nbasis();

          const double density_change_02 = max_density_change[i0 * size + i2];
          const double density_change_12 = max_density_change[i1 * size + i2];

          for (int i3 = i2; i3 != size; ++i3) {
            const unsigned int i23 = i2 * size + i3;
            if (i23 < i01) continue;
            int ijkl = plist->in_p4(i01, i23, i0, i1, i2, i3);
            if (ijkl == 0) continue;

            const double density_change_23 = max_density_change[i2 * size + i3] * 4.0;
            const double density_change_03 = max_density_change[i0 * size + i2];
            const double density_change_13 = max_density_change[i0 * size + i2];

            const bool eqli01i23 = (i01 == i23);

            const std::shared_ptr<const Shell>  b3 = basis[i3];
            const int b3offset = offset[i3];
            const int b3size = b3->nbasis();

            const double mulfactor = std::max(std::max(std::max(density_change_01, density_change_02),
                                             std::max(density_change_12, density_change_23)),
                                             std::max(density_change_03, density_change_13));
            const double integral_bound = mulfactor * schwarz_[i01] * schwarz_[i23];
            const bool skip_schwarz = integral_bound < schwarz_thresh_;
            if (skip_schwarz) continue;

            std::array<std::shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
#ifdef LIBINT_INTERFACE
            Libint eribatch(input);
#else
            ERIBatch eribatch(input, mulfactor);
#endif
            eribatch.compute();
            const double* eridata = eribatch.data();

            for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
              const int j0n = j0 * ndim_;

              for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                const unsigned int nj01 = (j0 << shift) + j1;
                const bool skipj0j1 = (j0 > j1);
                if (skipj0j1) {
                  eridata += b2size * b3size;
                  continue;
                }

                const bool eqlj0j1 = (j0 == j1);
                const double scal01 = (eqlj0j1 ? 0.5 : 1.0) * static_cast<double>(ijkl);
                const int j1n = j1 * ndim_;

                for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                  const int maxj1j2 = std::max(j1, j2);
                  const int minj1j2 = std::min(j1, j2);
                  const int minj1j2n = minj1j2 * ndim_;

                  const int maxj0j2 = std::max(j0, j2);
                  const int minj0j2 = std::min(j0, j2);
                  const int minj0j2n = minj0j2 * ndim_;
                  const int j2n = j2 * ndim_;

                  for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                    const bool skipj2j3 = (j2 > j3);
                    const unsigned int nj23 = (j2 << shift) + j3;
                    const bool skipj01j23 = (nj01 > nj23) && eqli01i23;

                    if (skipj2j3 || skipj01j23) continue;

                    const int maxj1j3 = std::max(j1, j3);
                    const int minj1j3 = std::min(j1, j3);

                    double intval = *eridata * scal01 * (j2 == j3 ? 0.5 : 1.0) * (nj01 == nj23 ? 0.25 : 0.5); // 1/2 in the Hamiltonian absorbed here
                    const double intval4 = 4.0 * intval;

                    data_[j0n + j1] += density_data[j2n + j3] * intval4;
                    data_[j2n + j3] += density_data[j0n + j1] * intval4;
                    data_[j0n + j3] -= density_data[j1n + j2] * intval;
                    data_[minj1j2n + maxj1j2] -= density_data[j0n + j3] * intval;
                    data_[minj0j2n + maxj0j2] -= density_data[j1n + j3] * intval;
                    data_[minj1j3 * ndim_ + maxj1j3] -= density_data[j0n + j2] * intval;
                  }
                }
              }
            }

          }
        }
      }
    }
    for (int i = 0; i != ndim_; ++i) data_[i*ndim_ + i] *= 2.0;
    fill_upper();

  //////////////////////////////////////////////////////////////////////////////////////////////////
  } else if (DF == 1) {

#ifndef NDEBUG
    std::cout << "    .. warning .. use a new Fock builder if possible (coeff_ required)" << std::endl;
#endif

    std::shared_ptr<const DFDist> df = geom_->df();

    // some constants
    const int naux = df->naux();
    assert(ndim_ == df->nbasis0());

    Timer pdebug(2);

    std::shared_ptr<Matrix> coeff = den_ex->copy();
    *coeff *= -1.0;
    int nocc = 0;
    {
      std::unique_ptr<double[]> vec(new double[ndim_]);
      coeff->diagonalize(vec.get());
      for (int i = 0; i != ndim_; ++i) {
        if (vec[i] < -1.0e-8) {
          ++nocc;
          dscal_(ndim_, std::sqrt(-vec[i]), coeff->data()+i*ndim_, 1);
        } else { break; }
      }
    }
    if (nocc == 0) return;
    pdebug.tick_print("Compute coeff (redundant)");

    std::shared_ptr<DFHalfDist> halfbj = df->compute_half_transform(coeff->slice(0,nocc));
    pdebug.tick_print("First index transform");

    std::shared_ptr<DFHalfDist> half = halfbj->apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->form_2index(half, -0.5);
    pdebug.tick_print("Exchange build");

    *this += *df->compute_Jop(density_);
    pdebug.tick_print("Coulomb build");
  }

}

template<int DF>
void Fock<DF>::fock_two_electron_part_with_coeff(const std::shared_ptr<const Matrix> ocoeff, const bool rhf, const double scale_exchange) {
  if (DF == 0) throw std::logic_error("Fock<DF>::fock_two_electron_part_with_coeff() is only for DF cases");

  Timer pdebug(2);

  std::shared_ptr<const DFDist> df = geom_->df();

  if (scale_exchange != 0.0) {
    std::shared_ptr<DFHalfDist> halfbj = df->compute_half_transform(ocoeff);
    pdebug.tick_print("First index transform");

    std::shared_ptr<DFHalfDist> half = halfbj->apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->form_2index(half, -1.0*scale_exchange);
    pdebug.tick_print("Exchange build");

    if (rhf) {
      auto coeff = std::make_shared<const Matrix>(*ocoeff->transpose()*2.0);
      *this += *df->compute_Jop(half, coeff, true);
    } else {
      *this += *df->compute_Jop(density_);
    }
  } else {
    *this += *df->compute_Jop(density_);
  }
  pdebug.tick_print("Coulomb build");
}

}


#endif
