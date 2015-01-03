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


#include <src/periodic/multipolebatch.h>
#include <src/integral/carsphlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/sortlist.h>
#include <src/integral/os/osintegral.h>

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

  // First, do (a|O_lm|0) using VRR
  const int worksize = num_multipoles_ * (amax_ + 1);
  complex<double>* const workx = stack_->get<complex<double>>(worksize);
  complex<double>* const worky = stack_->get<complex<double>>(worksize);
  complex<double>* const workz = stack_->get<complex<double>>(worksize);

  const int size_start = num_multipoles_ * asize_ * prim0_ * prim1_;
  complex<double>* intermediate_p = stack_->get<complex<double>>(size_start);

  const int nmul = num_multipoles_;

  for (int iprim = 0; iprim != prim0_ * prim1_; ++iprim) {

    vector<complex<double>*> current_data(nmul);
    for (int imul = 0; imul != nmul; ++imul)
      current_data[imul] = intermediate_p + iprim * asize_ + imul * asize_ * prim0_ * prim1_;

    for (int cnt = 0; cnt != nmul; ++cnt) {
      workx[cnt] = multipole_[iprim * prim0_ * prim1_ + cnt];
      worky[cnt] = 1.0;
      workz[cnt] = 1.0;
    }

    if (amax_ > 0) {
      int cnt = 0;
      for (int l = 0; l <= lmax_; ++l) {
        for (int m = 0; m <= 2*l; ++m, ++cnt) {
          workx[nmul + cnt] = - 2.0 * xb_[iprim] * AB_[0] * multipole_[iprim * prim0_ * prim1_ + cnt];
          worky[nmul + cnt] = - 2.0 * xb_[iprim] * AB_[1] * multipole_[iprim * prim0_ * prim1_ + cnt];
          workz[nmul + cnt] = - 2.0 * xb_[iprim] * AB_[2] * multipole_[iprim * prim0_ * prim1_ + cnt];
          if (abs(m - l) < l - 1) {
            const int i1 = (l - 1) * (l - 1) + m;
            workx[nmul + cnt] += 0.5 * (multipole_[iprim * prim0_ * prim1_ + i1] + multipole_[iprim * prim0_ * prim1_ + i1 - 2]);
            worky[nmul + cnt] += complex<double>(0.0, 0.5) * (multipole_[iprim * prim0_ * prim1_ + i1] + multipole_[iprim * prim0_ * prim1_ + i1 - 2]);
            workz[nmul + cnt] += multipole_[iprim * prim0_ * prim1_ + i1 - 1];
          }
        }
      }

      for (int a = 2; a <= amax_; ++a) {
        int cnt = 0;
        for (int l = 0; l <= lmax_; ++l) {
          for (int m = 0; m <= 2*l; ++m, ++cnt) {
            workx[a * nmul + cnt] = (a - 1.0) * workx[(a - 2) * nmul + cnt] - 2.0 * xb_[iprim] * AB_[0] * workx[(a - 1) * nmul + cnt];
            worky[a * nmul + cnt] = (a - 1.0) * worky[(a - 2) * nmul + cnt] - 2.0 * xb_[iprim] * AB_[1] * worky[(a - 1) * nmul + cnt];
            workz[a * nmul + cnt] = (a - 1.0) * workz[(a - 2) * nmul + cnt] - 2.0 * xb_[iprim] * AB_[2] * workz[(a - 1) * nmul + cnt];
            if (abs(m - l) < l - 1) {
              const int i1 = (l - 1) * (l - 1) + m;
              workx[a * nmul + cnt] += 0.5 * (workx[(a - 1) * nmul + i1] + multipole_[(a - 1) * nmul + i1 - 2]);
              worky[a * nmul + cnt] += complex<double>(0.0, 0.5) * (worky[(a - 1) * nmul + i1] + multipole_[(a - 1) * nmul + i1 - 2]);
              workz[a * nmul + cnt] += workz[(a - 1) * nmul + i1 - 1];
            }
          }
        }
      }
    }

    // get (a|O|0) from x, y, z components
    const int amax1 = amax_ + 1;
    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        for (int ix = 0; ix <= amax_ - iy - iz; ++ix) {

          for (int imul = 0; imul != nmul; ++imul) {
            const int iang = amapping_[ix + amax1 * (iy + amax1 * iz)];
            if (!swap01_) {
              current_data[imul][iang] = workx[ix * nmul + imul] * worky[iy * nmul + imul] * workz[iz * nmul + imul];
            } else {
              current_data[imul][iang] = -workx[ix * nmul + imul] * worky[iy * nmul + imul] * workz[iz * nmul + imul];
            }
          }

        }
      }
    }
  } // end loop over primitives

  const CSortList sort(spherical_);
  complex<double>* data_start = intermediate_p;
  complex<double>* data_final = data_;

  for (int imul = 0; imul != nmul; ++imul, data_start += imul * prim0_ * prim1_ * asize_,
                                           data_final += imul * cont0_ * cont1_ * asize_final_) {
    const int size_intermediate = asize_ * cont0_ * cont1_;
    complex<double>* intermediate_c = stack_->get<complex<double>>(size_intermediate);
    perform_contraction(asize_, data_start, prim0_, prim1_, intermediate_c,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    const int size_final_car = cont0_ * cont1_ * asize_intermediate_;
    complex<double>* const intermediate_fi = stack_->get<complex<double>>(size_final_car);

    // now get (a|O_lm|b) using HRR
    if (basisinfo_[1]->angular_number() != 0) {
      const int hrr_iprim = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr.hrrfunc_call(hrr_iprim, cont0_ * cont1_, intermediate_c, AB_, intermediate_fi);
    } else {
      copy_n(intermediate_c, size_final_car, intermediate_fi);
    }

    if (spherical_) {
      const int size_final_sph = cont0_ * cont1_ * asize_final_;
      complex<double>* const intermediate_i = stack_->get<complex<double>>(size_final_sph);
      const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      carsphlist.carsphfunc_call(carsph_index, cont0_ * cont1_, intermediate_fi, intermediate_i);

      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(sort_index, data_final, intermediate_i, cont1_, cont0_, 1, swap01_);
      stack_->release(size_final_sph, intermediate_i);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(sort_index, data_final, intermediate_fi, cont1_, cont0_, 1, swap01_);
    }

    stack_->release(size_final_car, intermediate_fi);
    stack_->release(size_intermediate, intermediate_c);
  } // end loop over multipoles


  // release resources
  stack_->release(size_start, intermediate_p);
  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}
