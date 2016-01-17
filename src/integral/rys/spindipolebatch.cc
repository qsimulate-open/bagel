//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spindipole.cc
// Copyright (C) 2015 Toru Shiozaki
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
#include <src/integral/rys/spin2rootlist.h>
#include <src/integral/rys/spindipolebatch.h>

using namespace std;
using namespace bagel;

const static HRRList hrr;
const static CarSphList carsphlist;

SpinDipoleBatch::SpinDipoleBatch(const array<shared_ptr<const Shell>,2>& _info, shared_ptr<const Atom> target, shared_ptr<StackMem> stack)
  : CoulombBatch_Base<double>(_info, make_shared<Molecule>(vector<shared_ptr<const Atom>>{target}, vector<shared_ptr<const Atom>>{}), 0, 2, stack), target_(target) {

  const double integral_thresh = PRIM_SCREEN_THRESH;
  this->allocate_arrays(primsize_);
  compute_ssss(integral_thresh);
  root_weight(primsize_);
}


void SpinDipoleBatch::compute() {

  double* const stack_save = stack_->template get<double>(size_block_);
  bkup_ = stack_save;

  const int worksize = rank_ * (amax1_+1);

  double* const workx = stack_->template get<double>(worksize);
  double* const worky = stack_->template get<double>(worksize);
  double* const workz = stack_->template get<double>(worksize);
  double* const worktx = stack_->template get<double>(worksize);
  double* const workty = stack_->template get<double>(worksize);
  double* const worktz = stack_->template get<double>(worksize);
  double* const worksx = stack_->template get<double>(worksize);
  double* const worksy = stack_->template get<double>(worksize);
  double* const worksz = stack_->template get<double>(worksize);

  const double Ax = basisinfo_[0]->position(0);
  const double Ay = basisinfo_[0]->position(1);
  const double Az = basisinfo_[0]->position(2);

  double r1x[20];
  double r1y[20];
  double r1z[20];
  double r2[20];

  fill_n(data_, size_alloc_, 0.0);

  const SortList sort(spherical1_);

  // perform VRR
  for (int xj = 0; xj != screening_size_; ++xj) {
    const int i = screening_[xj];
    const int offset_iprim = i * asize_;
    double* dataxx = data_  + offset_iprim;
    double* dataxy = dataxx + size_block_;
    double* dataxz = dataxy + size_block_;
    double* datayy = dataxz + size_block_;
    double* datayz = datayy + size_block_;
    double* datazz = datayz + size_block_;

    const double* croots = roots_ + i * rank_;
    const double* cweights = weights_ + i * rank_;
    for (int r = 0; r != rank_; ++r) {
      r1x[r] = P_[i*3  ] - Ax - (P_[i*3  ] - target_->position(0)) * croots[r];
      r1y[r] = P_[i*3+1] - Ay - (P_[i*3+1] - target_->position(1)) * croots[r];
      r1z[r] = P_[i*3+2] - Az - (P_[i*3+2] - target_->position(2)) * croots[r];
      r2[r] = (1.0 - croots[r]) * 0.5 / xp_[i];

      // absorb
      workx[r] = cweights[r] * coeff_[i] * xp_[i]*xp_[i] * 4.0 / (- target_->atom_charge());
      // not above that since this is reusing NAI code, Z is automatically multiplied and it has to be compensated here
      worky[r] = 1.0;
      workz[r] = 1.0;

      workx[rank_ + r] = r1x[r] * workx[r];
      worky[rank_ + r] = r1y[r] * worky[r];
      workz[rank_ + r] = r1z[r] * workz[r];
    }

    for (int j = 2; j <= amax1_; ++j) {
      const int offset = j * rank_;
      for (int r = 0; r != rank_; ++r) {
        workx[offset + r] = r1x[r] * workx[offset - rank_ + r] + (j - 1.0) * r2[r] * workx[offset - rank_ * 2 + r];
        worky[offset + r] = r1y[r] * worky[offset - rank_ + r] + (j - 1.0) * r2[r] * worky[offset - rank_ * 2 + r];
        workz[offset + r] = r1z[r] * workz[offset - rank_ + r] + (j - 1.0) * r2[r] * workz[offset - rank_ * 2 + r];
      }
    }

    // next compute tilde{I}_x, y, z in worktx,y,z (and remove 1-t from it - Eq. (26) of Shiozaki 2013)
    for (int j = 0; j <= amax1_; ++j) {
      const int offset = j * rank_;
      for (int r = 0; r != rank_; ++r) {
        worktx[offset + r] = (P_[i*3+0]-target_->position(0)) * workx[offset + r] + (j == 0 ? 0.0 : j*0.5/xp_[i] * workx[offset - rank_ + r]);
        workty[offset + r] = (P_[i*3+1]-target_->position(1)) * worky[offset + r] + (j == 0 ? 0.0 : j*0.5/xp_[i] * worky[offset - rank_ + r]);
        worktz[offset + r] = (P_[i*3+2]-target_->position(2)) * workz[offset + r] + (j == 0 ? 0.0 : j*0.5/xp_[i] * workz[offset - rank_ + r]);
      }
    }

    // then compute tilde{tilde{I}}_x,y,z
    for (int j = 0; j != amax1_; ++j) {
      const int offset = j * rank_;
      for (int r = 0; r != rank_; ++r) {
        worksx[offset + r] = worktx[offset + rank_ + r] + (Ax - target_->position(0)) * worktx[offset + r];
        worksy[offset + r] = workty[offset + rank_ + r] + (Ay - target_->position(1)) * workty[offset + r];
        worksz[offset + r] = worktz[offset + rank_ + r] + (Az - target_->position(2)) * worktz[offset + r];
      }
    }

    // assembly step
    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix) {
          const int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
          double xx = 0.0;
          double yy = 0.0;
          double zz = 0.0;
          for (int r = 0; r != rank_; ++r) {
            const double fac = 1.0 / 3.0 / (1.0-croots[r]);
            xx          += workz [iz * rank_ + r] * worky [iy * rank_ + r] * worksx[ix * rank_ + r] * fac;
            dataxy[pos] += workz [iz * rank_ + r] * workty[iy * rank_ + r] * worktx[ix * rank_ + r];
            dataxz[pos] += worktz[iz * rank_ + r] * worky [iy * rank_ + r] * worktx[ix * rank_ + r];
            yy          += workz [iz * rank_ + r] * worksy[iy * rank_ + r] * workx [ix * rank_ + r] * fac;
            datayz[pos] += worktz[iz * rank_ + r] * workty[iy * rank_ + r] * workx [ix * rank_ + r];
            zz          += worksz[iz * rank_ + r] * worky [iy * rank_ + r] * workx [ix * rank_ + r] * fac;
          }
          dataxx[pos] += 2.0*xx - yy - zz;
          datayy[pos] += 2.0*yy - zz - xx;
          datazz[pos] += 2.0*zz - xx - yy;
        }
      }
    }
  }

  // loop over cartesian blocks
  for (int i = 0; i != Nblocks(); ++i) {
    double* cdata = data_ + i * size_block_;
    fill_n(bkup_, size_block_, 0.0);

    // contract indices 01
    // data will be stored in bkup_: cont01{ xyz{ } }
    {
      const int m = asize_;
      this->perform_contraction(m, cdata, prim0size_, prim1size_, bkup_,
                          basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0size_,
                          basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1size_);
    }

    // HRR to indices 01
    // data will be stored in data_: cont01{ xyzab{ } }
    {
      if (basisinfo_[1]->angular_number() != 0) {
        const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
        hrr.hrrfunc_call(hrr_index, contsize_, bkup_, AB_, cdata);
      } else {
        copy_n(bkup_, size_block_, cdata);
      }
    }

    // Cartesian to spherical 01 if necesarry
    // data will be stored in bkup_
    if (spherical1_) {
      const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = contsize_;
      carsphlist.carsphfunc_call(carsphindex, nloops, cdata, bkup_);
    }

    // Sort cont01 and xyzab
    // data will be stored in data_: cont1b{ cont0a{ } }
    if (spherical1_) {
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(index, cdata, bkup_, cont1size_, cont0size_, 1, swap01_);
    } else {
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(index, bkup_, cdata, cont1size_, cont0size_, 1, swap01_);
      copy_n(bkup_, size_final_, cdata);
    }

  }

  stack_->release(worksize, worksz);
  stack_->release(worksize, worksy);
  stack_->release(worksize, worksx);
  stack_->release(worksize, worktz);
  stack_->release(worksize, workty);
  stack_->release(worksize, worktx);
  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
  stack_->release(size_block_, stack_save);
}


void SpinDipoleBatch::root_weight(const int ps) {
  spin2root__.root(rank_, T_, roots_, weights_, ps);
}


