//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: coulombbatch_energy.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/integral/sortlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/rys/coulombbatch_energy.h>

using namespace std;
using namespace bagel;

const static HRRList hrr;
const static CarSphList carsphlist;

void CoulombBatch_energy::compute() {

  double* const stack_save = stack_->template get<double>(size_block_);
  bkup_ = stack_save;

  const int worksize = rank_ * amax1_;

  double* const workx = stack_->template get<double>(worksize);
  double* const worky = stack_->template get<double>(worksize);
  double* const workz = stack_->template get<double>(worksize);

  const double Ax = basisinfo_[0]->position(0);
  const double Ay = basisinfo_[0]->position(1);
  const double Az = basisinfo_[0]->position(2);

  double r1x[20];
  double r1y[20];
  double r1z[20];
  double r2[20];

  fill_n(data_, size_block_, 0.0);

  const SortList sort(spherical1_);

  // perform VRR
  for (int xj = 0; xj != screening_size_; ++xj) {
    const int i = screening_[xj];
    const int ecp = (indexecp_.empty()) ? 1 : max_rterms_;
    const int iprim = i / (natom_ * ecp);
    const int iatom = (i / ecp) % natom_;
    const int offset_iprim = iprim * asize_;
    double* current_data = data_ + offset_iprim;

    const double* croots = roots_ + i * rank_;
    const double* cweights = weights_ + i * rank_;
    for (int r = 0; r != rank_; ++r) {
      double zeta = 1.0;
      if (!indexecp_.empty())
        zeta = mol_->atoms(iatom)->ecp_parameters()->shell_maxl_ecp()->ecp_exponents(indexecp_[xj]);
      const double sroot = scale_root(croots[r], xp_[i], zeta);
      const double sweight = scale_weight(cweights[r]);
      r1x[r] = P_[i * 3    ] - Ax - (P_[i * 3    ] - mol_->atoms(iatom)->position(0)) * sroot;
      r1y[r] = P_[i * 3 + 1] - Ay - (P_[i * 3 + 1] - mol_->atoms(iatom)->position(1)) * sroot;
      r1z[r] = P_[i * 3 + 2] - Az - (P_[i * 3 + 2] - mol_->atoms(iatom)->position(2)) * sroot;
      r2[r] = (1.0 - sroot) * 0.5 / xp_[i];

      workx[r] = sweight * coeff_[i];
      worky[r] = 1.0;
      workz[r] = 1.0;
    }

    if (amax_ > 0) {
      for (int r = 0; r != rank_; ++r) {
        workx[rank_ + r] = r1x[r] * workx[r];
        worky[rank_ + r] = r1y[r] * worky[r];
        workz[rank_ + r] = r1z[r] * workz[r];
      }
      for (int j = 2; j <= amax_; ++j) {
        const int offset = j * rank_;
        for (int r = 0; r != rank_; ++r) {
          workx[offset + r] = r1x[r] * workx[offset - rank_ + r] + (j - 1.0) * r2[r] * workx[offset - rank_ * 2 + r];
          worky[offset + r] = r1y[r] * worky[offset - rank_ + r] + (j - 1.0) * r2[r] * worky[offset - rank_ * 2 + r];
          workz[offset + r] = r1z[r] * workz[offset - rank_ + r] + (j - 1.0) * r2[r] * workz[offset - rank_ * 2 + r];
        }
      }
    }
    // assembly step
    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix) {
          const int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
          for (int r = 0; r != rank_; ++r)
          current_data[pos] += workz[iz * rank_ + r] * worky[iy * rank_ + r] * workx[ix * rank_ + r];
        }
      }
    }
  }

  // contract indices 01
  // data will be stored in bkup_: cont01{ xyz{ } }
  {
    const int m = asize_;
    this->perform_contraction(m, data_, prim0size_, prim1size_, bkup_,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0size_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1size_);
  }

  // HRR to indices 01
  // data will be stored in data_: cont01{ xyzab{ } }
  {
    if (basisinfo_[1]->angular_number() != 0) {
      const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr.hrrfunc_call(hrr_index, contsize_, bkup_, AB_, data_);
    } else {
      copy_n(bkup_, size_block_, data_);
    }
  }

  // Cartesian to spherical 01 if necesarry
  // data will be stored in bkup_
  if (spherical1_) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_;
    carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
  }

  // Sort cont01 and xyzab
  // data will be stored in data_: cont1b{ cont0a{ } }
  if (spherical1_) {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(index, data_, bkup_, cont1size_, cont0size_, 1, swap01_);
  } else {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(index, bkup_, data_, cont1size_, cont0size_, 1, swap01_);
    copy_n(bkup_, size_final_, data_);
  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
  stack_->release(size_block_, stack_save);
}

void CoulombBatch_energy::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (abs(T_[i]) < T_thresh__) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double erfsqt = inline_erf(sqrtt);
        weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    eriroot__.root(rank_, T_, roots_, weights_, ps);
  }
}


