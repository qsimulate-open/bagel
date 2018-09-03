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


void Orthogonal_Basis::add_shift(shared_ptr<const Orthogonal_Basis> t, const int istate, const double shift, const bool imag) {
  const double shift2 = shift * shift;

  for (int ist = 0; ist != nstates_; ++ist) {
    if (!sssr_ || ist == istate) {
      for (int iext = Excitations::aibj; iext != Excitations::total; ++iext) {
        const int pos = ist * Excitations::total + iext;
        shared_ptr<Tensor_<double>> ttensor = t->data(istate)->at(pos);
        const shared_ptr<Tensor_<double>> dtensor = denom_[istate]->at(pos);
        switch(iext) {
          case Excitations::aibj:
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0 : closed_) {
                    if (!dtensor->is_local(i0, i1, i2, i3)) continue;
                    unique_ptr<double[]> amplitude = ttensor->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> denom     = dtensor->get_block(i0, i1, i2, i3);
                    const size_t blocksize = ttensor->get_size(i0, i1, i2, i3);
                    if (imag) {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift2 / denom[j];
                      }
                    } else {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift;
                      }
                    }
                    data_[istate]->at(pos)->add_block(amplitude, i0, i1, i2, i3);
                  }
            break;
          case Excitations::arbs:
            for (auto& i3 : virt_)
              for (auto& i1 : virt_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i1, i3)) continue;
                  unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i1, i3);
                  unique_ptr<double[]> denom     = dtensor->get_block(i0o, i1, i3);
                  const size_t blocksize = ttensor->get_size(i0o, i1, i3);
                  if (imag) {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift2 / denom[j];
                    }
                  } else {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift;
                    }
                  }
                  data_[istate]->at(pos)->add_block(amplitude, i0o, i1, i3);
                }
            break;
          case Excitations::arbi:
            for (auto& i3 : virt_)
              for (auto& i2 : closed_)
                for (auto& i1 : virt_)
                  for (auto& i0o : interm_[iext]) {
                    if (!dtensor->is_local(i0o, i1, i2, i3)) continue;
                    unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i1, i2, i3);
                    unique_ptr<double[]> denom     = dtensor->get_block(i0o, i1, i2, i3);
                    const size_t blocksize = ttensor->get_size(i0o, i1, i2, i3);
                    if (imag) {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift2 / denom[j];
                      }
                    } else {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift;
                      }
                    }
                    data_[istate]->at(pos)->add_block(amplitude, i0o, i1, i2, i3);
                  }
            break;
          case Excitations::airj:
            for (auto& i2 : closed_)
              for (auto& i1 : virt_)
                for (auto& i0 : closed_)
                  for (auto& i0o : interm_[iext]) {
                    if (!dtensor->is_local(i0o, i0, i1, i2)) continue;
                    unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i0, i1, i2);
                    unique_ptr<double[]> denom     = dtensor->get_block(i0o, i0, i1, i2);
                    const size_t blocksize = ttensor->get_size(i0o, i0, i1, i2);
                    if (imag) {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift2 / denom[j];
                      }
                    } else {
                      for (size_t j = 0; j != blocksize; ++j) {
                        amplitude[j] *= shift;
                      }
                    }
                    data_[istate]->at(pos)->add_block(amplitude, i0o, i0, i1, i2);
                  }
            break;
          case Excitations::risj:
            for (auto& i2 : closed_)
              for (auto& i0 : closed_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i0, i2)) continue;
                  unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i0, i2);
                  unique_ptr<double[]> denom     = dtensor->get_block(i0o, i0, i2);
                  const size_t blocksize = ttensor->get_size(i0o, i0, i2);
                  if (imag) {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift2 / denom[j];
                    }
                  } else {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift;
                    }
                  }
                  data_[istate]->at(pos)->add_block(amplitude, i0o, i0, i2);
                }
            break;
          case Excitations::airs:
            for (auto& i1 : virt_)
              for (auto& i0 : closed_)
                for (auto& i0o : interm_[iext]) {
                  if (!dtensor->is_local(i0o, i0, i1)) continue;
                  unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i0, i1);
                  unique_ptr<double[]> denom     = dtensor->get_block(i0o, i0, i1);
                  const size_t blocksize = ttensor->get_size(i0o, i0, i1);
                  if (imag) {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift2 / denom[j];
                    }
                  } else {
                    for (size_t j = 0; j != blocksize; ++j) {
                      amplitude[j] *= shift;
                    }
                  }
                  data_[istate]->at(pos)->add_block(amplitude, i0o, i0, i1);
                }
            break;
          case Excitations::arst:
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i1)) continue;
                unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i1);
                unique_ptr<double[]> denom     = dtensor->get_block(i0o, i1);
                const size_t blocksize = ttensor->get_size(i0o, i1);
                if (imag) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift2 / denom[j];
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift;
                  }
                }
                data_[istate]->at(pos)->add_block(amplitude, i0o, i1);
              }
            break;
          case Excitations::rist:
            for (auto& i0 : closed_)
              for (auto& i0o : interm_[iext]) {
                if (!dtensor->is_local(i0o, i0)) continue;
                unique_ptr<double[]> amplitude = ttensor->get_block(i0o, i0);
                unique_ptr<double[]> denom     = dtensor->get_block(i0o, i0);
                const size_t blocksize = ttensor->get_size(i0o, i0);
                if (imag) {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift2 / denom[j];
                  }
                } else {
                  for (size_t j = 0; j != blocksize; ++j) {
                    amplitude[j] *= shift;
                  }
                }
                data_[istate]->at(pos)->add_block(amplitude, i0o, i0);
              }
            break;
        }
      }
    }
  }
}


tuple<shared_ptr<Matrix>,shared_ptr<Vec<double>>,shared_ptr<VecRDM<1>>,shared_ptr<VecRDM<2>>,shared_ptr<VecRDM<3>>,shared_ptr<VecRDM<3>>,vector<double>>
Orthogonal_Basis::make_d2_imag(shared_ptr<const Orthogonal_Basis> lambda, const double shift, const bool imag) const {
  auto dshift = make_shared<Matrix>(norb_, norb_);
  auto e0 = make_shared<Vec<double>>();
  auto e1 = make_shared<VecRDM<1>>();
  auto e2 = make_shared<VecRDM<2>>();
  auto e3 = make_shared<VecRDM<3>>();
  auto e4 = make_shared<VecRDM<3>>();
  vector<double> nimag;

  return tie(dshift, e0, e1, e2, e3, e4, nimag);
}

#endif
