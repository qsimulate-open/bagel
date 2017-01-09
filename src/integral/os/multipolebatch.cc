//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/integral/os/multipolebatch.h>
#include <src/integral/carsphlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/sortlist.h>

using namespace std;
using namespace bagel;

const static CHRRList hrr;
const static CCarSphList carsphlist;

MultipoleBatch::MultipoleBatch(const array<shared_ptr<const Shell>,2>& sh, const array<double, 3> centre,
                               const int lmax, shared_ptr<StackMem> stack)
 : MultipoleBatch_base(sh, centre, lmax, stack) {

 const double integral_thresh = PRIM_SCREEN_THRESH;
 allocate_arrays(prim0_ * prim1_);

 compute_ss(integral_thresh);
}

void MultipoleBatch::compute() {

  // First, do (ix, iy, iz; l; m) = (n|O_lm|0) using VRR

  complex<double>* intermediate_p = stack_->get<complex<double>>(size_alloc_);
  const complex<double> imag = swap01_ ? complex<double>(0.0, -1.0) : complex<double>(0.0, 1.0);

  const int nmul = num_multipoles_;
  vector<complex<double>*> bkup(nmul);
  for (int imul = 0; imul != nmul; ++imul)
    bkup[imul] = intermediate_p + imul * size_block_;

  for (int iprim = 0; iprim != prim0_ * prim1_; ++iprim) {

    vector<complex<double>> workz(amax1_ * nmul); // (0, 0, iz, l, m)

    // (0, 0, 0, l, m)
    for (int imul = 0; imul != nmul; ++imul)
      workz[imul] = multipole_[imul * prim0_ * prim1_ + iprim];

    int iang = 0;
    for (int iz = 0; iz <= amax_; ++iz) {
      // (0, 0, iz, l, m)
      if (iz > 0) {
        int imul = 0;
        for (int l = 0; l <= lmax_; ++l) {
          for (int m = 0; m <= 2 * l; ++m, ++imul) {
            assert(l * l + m == imul);
            workz[iz * nmul + imul] = -2.0 * xb_[iprim] * AB_[2] * workz[(iz-1)*nmul+imul];
            if (iz > 1)
              workz[iz * nmul + imul] += (iz-1.0) * workz[(iz-2)*nmul+imul];
            if (l > 0 && abs(m-l) < l) {
              const int i1 = (l-1)*(l-1) + m - 1;
              workz[iz * nmul + imul] += workz[(iz-1)*nmul+i1];
            }
            workz[iz * nmul + imul] *= 0.5 / xp_[iprim];
          }
        }
      }

      vector<complex<double>> worky((amax1_ - iz) * nmul);
      for (int imul = 0; imul != nmul; ++imul)
        worky[imul] = workz[iz * nmul + imul];

      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        // (0, iy, iz, l, m)
        for (int b = 1; b <= iy; ++b) {
          int imul = 0;
          for (int l = 0; l <= lmax_; ++l) {
            for (int m = 0; m <= 2 * l; ++m, ++imul) {
              worky[b * nmul + imul] = -2.0 * xb_[iprim] * AB_[1] * worky[(b-1)*nmul+imul];
              if (b > 1)
                worky[b * nmul + imul] += (b-1.0) * worky[(b-2)*nmul+imul];
              if (l > 0) {
                const int i1 = (l-1)*(l-1) + m - 1;
                if (abs(m - l - 1) < l)
                  worky[b * nmul + imul] -= imag * 0.5 * worky[(b-1)*nmul+i1-1];
                if (abs(m - l + 1) < l)
                  worky[b * nmul + imul] -= imag * 0.5 * worky[(b-1)*nmul+i1+1];
              }
              worky[b * nmul + imul] *= 0.5 / xp_[iprim];
            }
          }
        }

        vector<complex<double>> workx((amax1_ - iy - iz) * nmul);
        for (int imul = 0; imul != nmul; ++imul)
          workx[imul] = worky[iy * nmul + imul];

        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix, ++iang) {
          // (ix, iy, iz, l, m)
          for (int a = 1; a <= ix; ++a) {
            int imul = 0;
            for (int l = 0; l <= lmax_; ++l) {
              for (int m = 0; m <= 2 * l; ++m, ++imul) {
                workx[a * nmul + imul] = -2.0 * xb_[iprim] * AB_[0] * workx[(a-1)*nmul+imul];
                if (a > 1)
                  workx[a * nmul + imul] += (a-1.0) * workx[(a-2)*nmul+imul];
                if (l > 0) {
                  const int i1 = (l-1)*(l-1) + m - 1;
                  if (abs(m - l - 1) < l)
                    workx[a * nmul + imul] -= 0.5 * workx[(a-1)*nmul+i1-1];
                  if (abs(m - l + 1) < l)
                    workx[a * nmul + imul] += 0.5 * workx[(a-1)*nmul+i1+1];
                }
                workx[a * nmul + imul] *= 0.5 / xp_[iprim];
              }
            }
          }

          int imul = 0;
          for (int l = 0; l <= lmax_; ++l) {
            for (int m = 0; m <= 2 * l; ++m, ++imul) {
              complex<double>* current_data = bkup[imul] + iprim * asize_;
              const int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
              if (!swap01_) {
                current_data[pos] = workx[ix * nmul + imul];
              } else {
                //current_data[pos] = workx[ix * nmul + imul];
                current_data[pos] = conj(workx[ix * nmul + imul]);
              }
            }
          }
        } // ix
      } //iy
    } //iz
    assert(iang <= asize_);

  } // end loop over primitives - done (n|O|0)

  const CSortList sort(spherical_);
  complex<double>* data_start = intermediate_p;
  complex<double>* data_final = data_;
  fill_n(data_, size_alloc_, 0.0);


  for (int imul = 0; imul != nmul; ++imul, data_start += size_block_, data_final += size_block_) {

    complex<double>* const intermediate_c = stack_->get<complex<double>>(size_block_);
    perform_contraction(asize_, data_start, prim0_, prim1_, intermediate_c,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    complex<double>* const intermediate_fi = stack_->get<complex<double>>(size_block_);

    // now get (a|O_lm|b) using HRR
    if (basisinfo_[1]->angular_number() != 0) {
      const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr.hrrfunc_call(hrr_index, cont0_ * cont1_, intermediate_c, AB_, intermediate_fi);
    } else {
      copy_n(intermediate_c, size_block_, intermediate_fi);
    }

    if (spherical_) {
      complex<double>* const intermediate_i = stack_->get<complex<double>>(size_block_);
      const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      carsphlist.carsphfunc_call(carsph_index, cont0_ * cont1_, intermediate_fi, intermediate_i);

      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(sort_index, data_final, intermediate_i, cont1_, cont0_, 1, swap01_);
      stack_->release(size_block_, intermediate_i);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(sort_index, data_final, intermediate_fi, cont1_, cont0_, 1, swap01_);
    }

    stack_->release(size_block_, intermediate_fi);
    stack_->release(size_block_, intermediate_c);
  } // end loop over multipoles

  // release resources
  stack_->release(size_alloc_, intermediate_p);
}
