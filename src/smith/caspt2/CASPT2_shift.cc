//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_shift.cc
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

#include <src/smith/caspt2/CASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;


void CASPT2::CASPT2::add_shift(shared_ptr<VectorB> residual, shared_ptr<const VectorB> amplitude, const int istate) {
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = nact ? denom_->shalf_xx()->ndim()  * nvirt * nvirt : 0;
  const size_t size_arbi = nact ? denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt : 0;
  const size_t size_airj = nact ? denom_->shalf_h()->ndim()   * nclo * nvirt * nclo : 0;
  const size_t size_risj = nact ? denom_->shalf_hh()->ndim()  * nclo * nclo : 0;
  const size_t size_airs = nact ? denom_->shalf_xh()->ndim()  * nclo * nvirt : 0;
  const size_t size_arst = nact ? denom_->shalf_xxh()->ndim() * nvirt : 0;
  const size_t size_rist = nact ? denom_->shalf_xhh()->ndim() * nclo : 0;

  const double shift = info_->shift();
  const double shift2 = shift * shift;

  size_t ioffset = 0;

  for (int ist = 0; ist != nstates_; ++ist) {
    if (!info_->sssr() || ist == istate) {

      // a i b j case
      {
        const double e0loc = e0all_[ist] - e0_;
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                        const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = e0loc - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j3+nocc] + eig_[j1+nocc];
                        if (info_->shift_imag()) {
                          (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                        } else {
                          (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                        }
                      }
                }
              }
      }

      // a r b s case
      if (size_arbs) {
        ioffset += size_aibj;
        const size_t interm_size = denom_->shalf_xx()->ndim();
        for (auto& i3 : virt_)
          for (auto& i1 : virt_) {
            for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                for (int j02 = 0; j02 != interm_size; ++j02) {
                  const size_t jall = j02 + interm_size * (j1 + nvirt * j3);
                  const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j02) - e0all_[ist];
                  if (info_->shift_imag()) {
                    (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                  } else {
                    (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                  }
                }
              }
            }
          }
      }

      // a r b i case
      if (size_arbi) {
        ioffset += size_arbs;
        const size_t interm_size = denom_->shalf_x()->ndim();
        for (auto& i3 : virt_)
          for (auto& i2 : closed_)
            for (auto& i1 : virt_) {
              for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3)
                for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                    for (int j0 = 0; j0 != interm_size; ++j0) {
                      const size_t jall = j0 + interm_size * (j1 + nvirt * (j2 + nclo * j3));
                      const double denom = eig_[j1+nocc] + eig_[j3+nocc] - eig_[j2+ncore] + denom_->denom_x(j0) - e0all_[ist];
                      if (info_->shift_imag()) {
                        (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                      } else {
                        (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                      }
                    }
            }
      }

      // a i r j case
      if (size_airj) {
        ioffset += size_arbi;
        const size_t interm_size = denom_->shalf_h()->ndim();
        for (auto& i2 : closed_)
          for (auto& i1 : virt_)
            for (auto& i0 : closed_) {
              for (int j3 = 0; j3 != interm_size; ++j3)
                for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                    for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                      const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                      const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j3) - e0all_[ist];
                      if (info_->shift_imag()) {
                        (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                      } else {
                        (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                      }
                    }
            }
      }

      // r i s j case
      if (size_risj) {
        ioffset += size_airj;
        const size_t interm_size = denom_->shalf_hh()->ndim();
        for (auto& i2 : closed_)
          for (auto& i0 : closed_) {
            for (int j13 = 0; j13 != interm_size; ++j13)
              for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                  const size_t jall = j0 + nclo * (j2 + nclo * j13);
                  const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_hh(j13) - e0all_[ist];
                  if (info_->shift_imag()) {
                    (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                  } else {
                    (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                  }
                }
          }
      }

      // a i r s & a r s i case
      // TODO implement complex case
      if (size_airs) {
        ioffset += size_risj;
        const size_t interm_size = denom_->shalf_xh()->ndim();
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            for (int j23 = 0; j23 != interm_size; ++j23)
              for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                  const size_t jall = j0 + nclo * (j1 + nvirt * j23);
                  const double denom = eig_[j1+nocc] - eig_[j0+ncore] + denom_->denom_xh(j23) - e0all_[ist];
                  if (info_->shift_imag()) {
                    (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                  } else {
                    (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                  }
               }
          }
      }

      // a r s t case
      if (size_arst) {
        ioffset += size_airs;
        const size_t interm_size = denom_->shalf_xxh()->ndim();
        for (auto& i1 : virt_) {
          for (int j023 = 0; j023 != interm_size; ++j023)
            for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
              const size_t jall = j1 + nvirt * j023;
              const double denom = eig_[j1+nocc] + denom_->denom_xxh(j023) - e0all_[ist];
              if (info_->shift_imag()) {
                (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
              } else {
                (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
              }
            }
        }
      }

      // r i s t case
      if (size_rist) {
        ioffset += size_arst;
        const size_t interm_size = denom_->shalf_xhh()->ndim();
        for (auto& i2 : closed_) {
          for (int j013 = 0; j013 != interm_size; ++j013)
            for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2) {
              const size_t jall = j2 + nclo * j013;
              const double denom = - eig_[j2+ncore] + denom_->denom_xhh(j013) - e0all_[ist];
              if (info_->shift_imag()) {
                (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
              } else {
                (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
              }
            }
        }
        ioffset += size_rist;
      }
    }
  }
}

#endif
