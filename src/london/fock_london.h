//
// BAGEL - Parallel electron correlation program.
// Filename: fock_london.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __BAGEL_SRC_LONDON_FOCK_LONDON_H
#define __BAGEL_SRC_LONDON_FOCK_LONDON_H

#include <src/df/df_london.h>
#include <src/integral/libint/libint.h>
#include <src/integral/comprys/complexeribatch.h>
#include <src/london/scf_base_london.h>
#include <src/london/fock_base_london.h>

namespace bagel {

template<int DF>
class Fock_London : public Fock_base_London {
  protected:
    void fock_two_electron_part(std::shared_ptr<const ZMatrix> den = nullptr);
    void fock_two_electron_part_with_coeff(const std::shared_ptr<const ZMatrix> coeff, const bool rhf, const double scale_ex);

    // when DF gradients are requested
    bool store_half_;
    std::shared_ptr<DFHalfDist_London> half_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Fock_base_London>(*this) & store_half_;
    }

  public:
    Fock_London() { }
    // Fock operator for DF cases
    template<int DF1 = DF, class = typename std::enable_if<DF1==1>::type>
    Fock_London(const std::shared_ptr<const Geometry_London> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c,
         const std::shared_ptr<const ZMatrix> ocoeff, const bool store = false, const bool rhf = false, const double scale_ex = 1.0)
     : Fock_base_London(a,b,c), store_half_(store) {
      fock_two_electron_part_with_coeff(ocoeff, rhf, scale_ex);
      fock_one_electron_part();
    }

    // Fock operator
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock_London(const std::shared_ptr<const Geometry_London> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c, const std::vector<double>& d) : Fock_London(a,b,c,c,d) {}

    // Fock operator with a different density matrix for exchange
    template<int DF1 = DF, class = typename std::enable_if<DF1==1 or DF1==0>::type>
    Fock_London(const std::shared_ptr<const Geometry_London> a, const std::shared_ptr<const ZMatrix> b, const std::shared_ptr<const ZMatrix> c, std::shared_ptr<const ZMatrix> ex,
         const std::vector<double>& d)
     : Fock_base_London(a,b,c,d), store_half_(false) {
      fock_two_electron_part(ex);
      fock_one_electron_part();
    }

    std::shared_ptr<DFHalfDist_London> half() const { return half_; }
};


template<int DF>
void Fock_London<DF>::fock_two_electron_part(std::shared_ptr<const ZMatrix> den_ex) {

  const std::vector<std::shared_ptr<const Atom>> atoms = cgeom_->atoms();
  std::vector<std::shared_ptr<const Shell>> basis;
  std::vector<int> offset;
  int cnt = 0;
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++cnt) {
    const std::vector<std::shared_ptr<const Shell>> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());
    const std::vector<int> tmpoff = cgeom_->offset(cnt);
    offset.insert(offset.end(), tmpoff.begin(), tmpoff.end());
  }

  //const int shift = sizeof(int) * 4;
  const int size = basis.size();

  // first make max_density_change std::vector for each batch pair.
  const std::complex<double>* density_data = density_->data();

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
  if (DF == 0) {
    //////////////// ONLY FOR REFERENCES. //////////////////
#if 1
    std::shared_ptr<Petite> plist = cgeom_->plist();;

    for (int i0 = 0; i0 != size; ++i0) {
      //if (!plist->in_p1(i0)) continue;

      const std::shared_ptr<const Shell>  b0 = basis[i0];
      const int b0offset = offset[i0];
      const int b0size = b0->nbasis();
      for (int i1 = 0; i1 != size; ++i1) {
      //for (int i1 = i0; i1 != size; ++i1) {
        const unsigned int i01 = i0 *size + i1;
        //if (!plist->in_p2(i01)) continue;

        const std::shared_ptr<const Shell>  b1 = basis[i1];
        const int b1offset = offset[i1];
        const int b1size = b1->nbasis();

        const double density_change_01 = max_density_change[i01] * 4.0;

        for (int i2 = 0; i2 != size; ++i2) {
        //for (int i2 = i0; i2 != size; ++i2) {
          const std::shared_ptr<const Shell>  b2 = basis[i2];
          const int b2offset = offset[i2];
          const int b2size = b2->nbasis();

          const double density_change_02 = max_density_change[i0 * size + i2];
          const double density_change_12 = max_density_change[i1 * size + i2];

          for (int i3 = 0; i3 != size; ++i3) {
          //for (int i3 = i2; i3 != size; ++i3) {
            const unsigned int i23 = i2 * size + i3;
            //if (i23 < i01) continue;
            //int ijkl = plist->in_p4(i01, i23, i0, i1, i2, i3);
            //if (ijkl == 0) continue;

            const double density_change_23 = max_density_change[i2 * size + i3] * 4.0;
            const double density_change_03 = max_density_change[i0 * size + i2];
            const double density_change_13 = max_density_change[i0 * size + i2];

            //const bool eqli01i23 = (i01 == i23);

            const std::shared_ptr<const Shell>  b3 = basis[i3];
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

            std::array<std::shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
            ComplexERIBatch eribatch(input, mulfactor);
            eribatch.compute();
            const std::complex<double>* eridata = eribatch.data();
            for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
              const int j0n = j0 * ndim_;

              for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
                //const unsigned int nj01 = (j0 << shift) + j1;
                //const bool skipj0j1 = (j0 > j1);
                //if (skipj0j1) {
                //  eridata += b2size * b3size;
                //  continue;
                //}

                //const bool eqlj0j1 = (j0 == j1);
                //const double scal01 = (eqlj0j1 ? 0.5 : 1.0) * static_cast<double>(ijkl);
                const int j1n = j1 * ndim_;

                for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                  //const int maxj1j2 = std::max(j1, j2);
                  //const int minj1j2 = std::min(j1, j2);
                  //const int minj1j2n = minj1j2 * ndim_;

                  //const int maxj0j2 = std::max(j0, j2);
                  //const int minj0j2 = std::min(j0, j2);
                  //const int minj0j2n = minj0j2 * ndim_;
                  const int j2n = j2 * ndim_;

                  for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                    //const bool skipj2j3 = (j2 > j3);
                    //const unsigned int nj23 = (j2 << shift) + j3;
                    //const bool skipj01j23 = (nj01 > nj23) && eqli01i23;

                    //if (skipj2j3 || skipj01j23) continue;

                    //const int maxj1j3 = std::max(j1, j3);
                    //const int minj1j3 = std::min(j1, j3);

                    //std::complex<double> intval = *eridata * scal01 * (j2 == j3 ? 0.5 : 1.0) * (nj01 == nj23 ? 0.25 : 0.5); // 1/2 in the Hamiltonian absorbed here
                    //const std::complex<double> intval4 = 4.0 * intval;

                    const int j3n = j3 * ndim_;
                    std::complex<double> intval = *eridata * 0.5; // 1/2 in the Hamiltonian absorbed here

                    if (j0 + j1 <  j2 + j3) continue;
                    if (j0 + j2 <  j1 + j3) continue;
                    if (j0 + j1 == j2 + j3) intval *= 0.5;
                    if (j0 + j2 == j1 + j3) intval *= 0.5;
                    const std::complex<double> intval2 = intval * 2.0;

                    data_[j0n + j1] += density_data[j3n + j2] * intval2; // Coulomb  (ab|cd)
                    data_[j0n + j3] -= density_data[j1n + j2] * intval;  // Exchange (ad|cb)

                    data_[j2n + j3] += density_data[j1n + j0] * intval2; // Coulomb  (cd|ab)
                    data_[j2n + j1] -= density_data[j3n + j0] * intval;  // Exchange (cb|ad)

                    data_[j1n + j0] += density_data[j2n + j3] * std::conj(intval2); // Coulomb  (ba|dc)
                    data_[j3n + j0] -= density_data[j2n + j1] * std::conj(intval);  // Exchange (da|bc)

                    data_[j3n + j2] += density_data[j0n + j1] * std::conj(intval2); // Coulomb  (dc|ba)
                    data_[j1n + j2] -= density_data[j0n + j3] * std::conj(intval);  // Exchange (bc|dz)

                    /*
                    data_[j0n + j1] += density_data[j2n + j3] * intval4;
                    data_[j2n + j3] += density_data[j0n + j1] * intval4;
                    data_[j0n + j3] -= density_data[j1n + j2] * intval;
                    data_[minj1j2n + maxj1j2] -= density_data[j0n + j3] * intval;
                    data_[minj0j2n + maxj0j2] -= density_data[j1n + j3] * intval;
                    data_[minj1j3 * ndim_ + maxj1j3] -= density_data[j0n + j2] * intval;
                    */
                  }
                }
              }
            }

          }
        }
      }
    }
    //for (int i = 0; i != ndim_; ++i) data_[i*ndim_ + i] *= 2.0;

#endif
  //////////////////////////////////////////////////////////////////////////////////////////////////
  } else if (DF == 1) {

#ifndef NDEBUG
    std::cout << "    .. warning .. use a new Fock builder if possible (coeff_ required)" << std::endl;
#endif

    std::shared_ptr<const DFDist_London> df = cgeom_->df();

    // some constants
    assert(ndim_ == df->nbasis0());

    Timer pdebug(3);

    std::shared_ptr<ZMatrix> coeff = den_ex->copy();
    *coeff *= -1.0;
    int nocc = 0;
    {
      std::unique_ptr<double[]> vec(new double[ndim_]);
      coeff->diagonalize(vec.get());
      for (int i = 0; i != ndim_; ++i) {
        if (vec[i] < -1.0e-8) {
          ++nocc;
          const double fac = std::sqrt(-vec[i]);
          std::for_each(coeff->element_ptr(0,i), coeff->element_ptr(0,i+1), [&fac](std::complex<double>& i) { i *= fac; });
        } else { break; }
      }
    }
    if (nocc == 0) return;
    pdebug.tick_print("Compute coeff (redundant)");

    std::shared_ptr<DFHalfDist_London> halfbj = df->compute_half_transform(coeff->slice(0,nocc));
    pdebug.tick_print("First index transform");

    std::shared_ptr<DFHalfDist_London> half = halfbj->apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->form_2index(half, -0.5);
    pdebug.tick_print("Exchange build");

    *this += *df->compute_Jop(density_);
    pdebug.tick_print("Coulomb build");
  }

}

template<int DF>
void Fock_London<DF>::fock_two_electron_part_with_coeff(const std::shared_ptr<const ZMatrix> ocoeff, const bool rhf, const double scale_exchange) {
  if (DF == 0) throw std::logic_error("Fock_London<DF>::fock_two_electron_part_with_coeff() is only for DF cases");

#if 1
  Timer pdebug(3);

  std::shared_ptr<const DFDist_London> df = cgeom_->df();

  if (scale_exchange != 0.0) {
    std::shared_ptr<DFHalfDist_London> halfbj = df->compute_half_transform(ocoeff);
    pdebug.tick_print("First index transform");

    std::shared_ptr<DFHalfDist_London> half = halfbj->apply_J();
    pdebug.tick_print("Metric multiply");

    *this += *half->form_2index(half, -1.0*scale_exchange);
    pdebug.tick_print("Exchange build");

    if (rhf) {
      auto coeff = std::make_shared<const ZMatrix>(*ocoeff->transpose()*2.0);
      *this += *df->compute_Jop(half, coeff, true);
    } else {
      *this += *df->compute_Jop(density_);
      assert(0);
    }
    // when gradient is requested..
    if (store_half_)
      half_ = half;
  } else {
    *this += *df->compute_Jop(density_);
    assert(0);
  }
  pdebug.tick_print("Coulomb build");
#endif

}

}

extern template class bagel::Fock_London<0>;
extern template class bagel::Fock_London<1>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Fock_London<0>)
BOOST_CLASS_EXPORT_KEY(bagel::Fock_London<1>)

#endif
