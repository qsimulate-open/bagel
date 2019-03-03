//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: orthogonal_shift.cc
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <numeric>
#include <src/smith/moint.h>
#include <src/smith/orthogonal.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


void Orthogonal_Basis::add_shift(shared_ptr<const Orthogonal_Basis> t, const int istate) {
  const double shift2 = shift_ * shift_;

  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    const shared_ptr<Tensor_<double>> ttensor = t->data(istate)->at(iext);
    const int dataindex = sssr_ ? iext + istate * Excitations::total : iext;
    switch(iext) {
      case Excitations::arbs:
        for (auto& i3 : virt_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[dataindex]) {
              if (!ttensor->is_local(i0o, i1, i3)) continue;
              unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i1, i3);
              const size_t blocksize = ttensor->get_size(i0o, i1, i3);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift2;
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift_;
                }
              }
              data_[istate]->at(iext)->add_block(amplitude, i0o, i1, i3);
            }
        break;
      case Excitations::arbi:
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[dataindex]) {
                if (!ttensor->is_local(i0o, i1, i2, i3)) continue;
                unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i1, i2, i3);
                const size_t blocksize = ttensor->get_size(i0o, i1, i2, i3);
                if (imag_) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift2;
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift_;
                  }
                }
                data_[istate]->at(iext)->add_block(amplitude, i0o, i1, i2, i3);
              }
        break;
      case Excitations::airj:
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[dataindex]) {
                if (!ttensor->is_local(i0o, i0, i1, i2)) continue;
                unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i0, i1, i2);
                const size_t blocksize = ttensor->get_size(i0o, i0, i1, i2);
                if (imag_) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift2;
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift_;
                  }
                }
                data_[istate]->at(iext)->add_block(amplitude, i0o, i0, i1, i2);
              }
        break;
      case Excitations::risj:
        for (auto& i2 : closed_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[dataindex]) {
              if (!ttensor->is_local(i0o, i0, i2)) continue;
              unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i0, i2);
              const size_t blocksize = ttensor->get_size(i0o, i0, i2);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift2;
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift_;
                }
              }
              data_[istate]->at(iext)->add_block(amplitude, i0o, i0, i2);
            }
        break;
      case Excitations::airs:
        for (auto& i1 : virt_)
          for (auto& i0 : closed_)
            for (auto& i0o : interm_[dataindex]) {
              if (!ttensor->is_local(i0o, i0, i1)) continue;
              unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i0, i1);
              const size_t blocksize = ttensor->get_size(i0o, i0, i1);
              if (imag_) {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift2;
                }
              } else {
                for (size_t j = 0; j != blocksize; ++j) {
                  amplitude[j] *= shift_;
                }
              }
              data_[istate]->at(iext)->add_block(amplitude, i0o, i0, i1);
            }
        break;
      case Excitations::arst:
        for (auto& i1 : virt_)
          for (auto& i0o : interm_[dataindex]) {
            if (!ttensor->is_local(i0o, i1)) continue;
            unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i1);
            const size_t blocksize = ttensor->get_size(i0o, i1);
            if (imag_) {
              for (size_t j = 0; j != blocksize; ++j) {
                amplitude[j] *= shift2;
              }
            } else {
              for (size_t j = 0; j != blocksize; ++j) {
                amplitude[j] *= shift_;
              }
            }
            data_[istate]->at(iext)->add_block(amplitude, i0o, i1);
          }
        break;
      case Excitations::rist:
        for (auto& i0 : closed_)
          for (auto& i0o : interm_[dataindex]) {
            if (!ttensor->is_local(i0o, i0)) continue;
            unique_ptr<double[]> amplitude   = ttensor->get_block(i0o, i0);
            const size_t blocksize = ttensor->get_size(i0o, i0);
            if (imag_) {
              for (size_t j = 0; j != blocksize; ++j) {
                amplitude[j] *= shift2;
              }
            } else {
              for (size_t j = 0; j != blocksize; ++j) {
                amplitude[j] *= shift_;
              }
            }
            data_[istate]->at(iext)->add_block(amplitude, i0o, i0);
          }
        break;
      case Excitations::aibj:
        for (int ist = 0; ist != nstates_; ++ist) {
          if (!sssr_ || ist == istate) {
            const int pos = iext + (sssr_ ? 0 : ist);
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!ttensor->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> amplitude   = t->data(istate)->at(pos)->get_block(i0, i1, i2, i3);
                    const size_t blocksize = ttensor->get_size(i0, i1, i2, i3);
                    if (imag_) {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift2;
                      }
                    } else {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift_;
                      }
                    }
                    data_[istate]->at(pos)->add_block(amplitude, i0, i1, i2, i3);
                  }
          }
        }
        break;
    }
  }
  mpi__->barrier();
}

#endif
