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
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

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
      {
        ioffset += size_aibj;
        for (auto& i2 : active_)
          for (auto& i0 : active_) {
            const size_t interm_size = denom_->shalf_xx()->ndim();
            for (auto& i3 : virt_)
              for (auto& i1 : virt_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3) {
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                    for (int j02 = 0; j02 != interm_size; ++j02) {
                      const size_t jall = j02 + interm_size * (j1 + nvirt * j3);
                      const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j02) - e0_;
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
      }

      // a r b i case
      {
        ioffset += size_arbs;
        for (auto& i0 : active_) {
          const size_t interm_size = denom_->shalf_x()->ndim();
          for (auto& i3 : virt_)
            for (auto& i2 : closed_)
              for (auto& i1 : virt_) {
                for (int j3 = i3.offset()-nocc; j3 != i3.offset()+i3.size()-nocc; ++j3)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = 0; j0 != interm_size; ++j0) {
                        const size_t jall = j0 + interm_size * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = eig_[j1+nocc] + eig_[j3+nocc] - eig_[j2+ncore] + denom_->denom_x(j0) - e0_;
                        if (info_->shift_imag()) {
                          (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                        } else {
                          (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                        }
                      }
              }
        }
      }

      // a i r j case
      {
        ioffset += size_arbi;
        for (auto& i3 : active_) {
          const size_t interm_size = denom_->shalf_h()->ndim();
          for (auto& i2 : closed_)
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j3 = 0; j3 != interm_size; ++j3)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                      for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                        const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3));
                        const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j3) - e0_;
                        if (info_->shift_imag()) {
                          (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                        } else {
                          (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                        }
                      }
              }
        }
      }

      // r i s j case
      {
        ioffset += size_airj;
        for (auto& i3 : active_)
          for (auto& i1 : active_) {
            const size_t interm_size = denom_->shalf_hh()->ndim();
            for (auto& i2 : closed_)
              for (auto& i0 : closed_) {
                for (int j13 = 0; j13 != interm_size; ++j13)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2)
                    for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                      const size_t jall = j0 + nclo * (j2 + nclo * j13);
                      const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_hh(j13) - e0_;
                      if (info_->shift_imag()) {
                        (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                      } else {
                        (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                      }
                    }
              }
          }
      }

      // a i r s & a r s i case
      // TODO implement complex case
      {
        ioffset += size_risj;
        for (auto& i3 : active_)
          for (auto& i2 : active_) {
            const size_t interm_size = denom_->shalf_xh()->ndim();
            for (auto& i1 : virt_)
              for (auto& i0 : closed_) {
                for (int j23 = 0; j23 != interm_size; ++j23)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1)
                    for (int j0 = i0.offset()-ncore; j0 != i0.offset()+i0.size()-ncore; ++j0) {
                      const size_t jall = j0 + nclo * (j1 + nvirt * j23);
                      const double denom = eig_[j1+nocc] - eig_[j0+ncore] + denom_->denom_xh(j23) - e0_;
                      if (info_->shift_imag()) {
                        (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                      } else {
                        (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                      }
                   }
              }
          }
      }

      // a r s t case
      {
        ioffset += size_airs;
        for (auto& i3 : active_)
          for (auto& i2 : active_)
            for (auto& i0 : active_) {
              const size_t interm_size = denom_->shalf_xxh()->ndim();
              for (auto& i1 : virt_) {
                for (int j023 = 0; j023 != interm_size; ++j023)
                  for (int j1 = i1.offset()-nocc; j1 != i1.offset()+i1.size()-nocc; ++j1) {
                    const size_t jall = j1 + nvirt * j023;
                    const double denom = eig_[j1+nocc] + denom_->denom_xxh(j023) - e0_;
                    if (info_->shift_imag()) {
                      (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                    } else {
                      (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                    }
                  }
              }
            }
      }

      // r i s t case
      {
        ioffset += size_arst;
        for (auto& i3 : active_)
          for (auto& i1 : active_)
            for (auto& i0 : active_) {
              const size_t interm_size = denom_->shalf_xhh()->ndim();
              for (auto& i2 : closed_) {
                for (int j013 = 0; j013 != interm_size; ++j013)
                  for (int j2 = i2.offset()-ncore; j2 != i2.offset()+i2.size()-ncore; ++j2) {
                    const size_t jall = j2 + nclo * j013;
                    const double denom = - eig_[j2+ncore] + denom_->denom_xhh(j013) - e0_;
                    if (info_->shift_imag()) {
                      (*residual)[ioffset + jall] += shift2 * (*amplitude)[ioffset + jall] / denom;
                    } else {
                      (*residual)[ioffset + jall] += shift * (*amplitude)[ioffset + jall];
                    }
                  }
              }
            }
        ioffset += size_rist;
      }

    }
  }
}


shared_ptr<Matrix> CASPT2::CASPT2::make_d2_imag(vector<shared_ptr<VectorB>> lambda, vector<shared_ptr<VectorB>> amplitude) const {
  shared_ptr<Matrix> dshift = den2_->clone();
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;
  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;
  const double shift2 = info_->shift() * info_->shift();
  for (int istate = 0; istate != nstates_; ++istate) { // state of T
    // temporary. will replace (RDM calculation redundant)
    shared_ptr<const RDM<1>> rdm1;
    shared_ptr<const RDM<2>> rdm2;
    shared_ptr<const RDM<3>> rdm3;
    shared_ptr<const RDM<4>> rdm4;

    // considering only SS-SR case
    tie(rdm1, rdm2) = info_->rdm12(istate, istate);
    tie(rdm3, rdm4) = info_->rdm34(istate, istate);
    shared_ptr<const VectorB> l = lambda[istate];
    shared_ptr<const VectorB> t = amplitude[istate];
    // a i b j
    size_t ioffset = 0;
    {
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j2 = 0; j2 != nclo; ++j2)
          for (size_t j1 = 0; j1 != nvirt; ++j1)
            for (size_t j0 = 0; j0 != nclo; ++j0) {
              const size_t j0i = j0;
              const size_t j1i = j1 + nocc - ncore;
              const size_t j2i = j2;
              const size_t j3i = j3 + nocc - ncore;
              const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
              const size_t jall2 = j0 + nclo * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
              const double lcovar = ((*l)[jall] * 8.0 - (*l)[jall2] * 4.0);
              const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j1+nocc] + eig_[j3+nocc];
              const double z = lcovar * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j0i, j0i) += z;
              dshift->element(j1i, j1i) -= z;
              dshift->element(j2i, j2i) += z;
              dshift->element(j3i, j3i) -= z;
            }
      ioffset += size_aibj;
    }

    // a r b s
    {
      const size_t interm_size = denom_->shalf_xx()->ndim();
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j1 = 0; j1 != nvirt; ++j1) {
            const size_t j1i = j1 + nocc - ncore;
            const size_t j3i = j3 + nocc - ncore;
            for (size_t j02 = 0; j02 != interm_size; ++j02) {
              const size_t jall = j02 + interm_size * (j1 + nvirt * j3) + ioffset;
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j02) - e0_;
              const double z = (*l)[jall] * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j1i, j1i) -= z;
              dshift->element(j3i, j3i) -= z;
              for (size_t j2 = 0; j2 != nact; ++j2) {
                const size_t j2i = j2 + nclo;
                for (size_t j0 = 0; j0 != nact; ++j0) {
                  const size_t j0i = j0 + nclo;
                  for (size_t j4 = 0; j4 != nact; ++j4)
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      const double VtuO = denom_->shalf_xx()->element(j02, j4+j5*nact + istate*nact*nact);
                      for (size_t j6 = 0; j6 != nact; ++j6)
                        for (size_t j7 = 0; j7 != nact; ++j7) {
                          const double VvwO = denom_->shalf_xx()->element(j02, j6+j7*nact + istate*nact*nact);
                          dshift->element(j0i, j2i) -= z * rdm3->element(j4, j5, j6, j7, j0, j2) * VtuO * VvwO;
                        }
                    }
                  dshift->element(j0i, j2i) += z * rdm1->element(j0, j2);
                }
              }
            }
        }
      ioffset += size_arbs;
    }

    // a r b i
    {
      const size_t interm_size = denom_->shalf_x()->ndim();
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j2 = 0; j2 != nclo; ++j2)
          for (size_t j1 = 0; j1 != nvirt; ++j1) {
            const size_t j3i = j3 + nocc - ncore;
            const size_t j2i = j2;
            const size_t j1i = j1 + nocc - ncore;
            for (size_t j0o = 0; j0o != interm_size; ++j0o) {
              const size_t jall = j0o + nclo * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
              const size_t jall2 = j0o + nclo * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j0o) - e0_;
              const double z = lcovar * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j1i, j1i) -= z;
              dshift->element(j3i, j3i) -= z;
              dshift->element(j2i, j2i) += z;
              for (size_t j0 = 0; j0 != nact; ++j0) {
                const size_t j0i = j0 + nclo;
                for (size_t j6 = 0; j6 != nact; ++j6) {
                  const size_t j6i = j6 + nclo;
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    const double VtO = denom_->shalf_x()->element(j0o, j4 + istate*nact);
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      const double VuO = denom_->shalf_x()->element(j0o, j5 + istate*nact);
                      dshift->element(j0i, j6i) -= z * rdm2->element(j4, j5, j0, j6) * VtO * VuO;
                    }
                  }
                  dshift->element(j0i, j6i) += z * rdm1->element(j0, j6);
                }
              }
            }
          }
      ioffset += size_arbi;
    }
    
    // a i r j
    {
      const size_t interm_size = denom_->shalf_h()->ndim();
      for (size_t j2 = 0; j2 != nclo; ++j2)
        for (size_t j1 = 0; j1 != nvirt; ++j1)
          for (size_t j0 = 0; j0 != nclo; ++j0) {
            const size_t j2i = j2;
            const size_t j1i = j1 + nocc - ncore;
            const size_t j0i = j0;
            for (size_t j3o = 0; j3o != interm_size; ++j3o) {
              const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3o)) + ioffset;
              const size_t jall2 = j2 + nclo * (j1 + nvirt * (j0 + nclo * j3o)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j3o) - e0_;
              const double z = lcovar * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j1i, j1i) -= z;
              dshift->element(j0i, j0i) += z;
              dshift->element(j2i, j2i) += z;
              for (size_t j3 = 0; j3 != nact; ++j3) {
                const size_t j3i = j3 + nclo;
                for (size_t j6 = 0; j6 != nact; ++j6) {
                  const size_t j6i = j6 + nclo;
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    const double VtO = denom_->shalf_h()->element(j3o, j4 + istate*nact);
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      const double VuO = denom_->shalf_h()->element(j3o, j5 + istate*nact);
                      dshift->element(j3i, j6i) += z * rdm2->element(j4, j5, j3, j6) * VtO * VuO;
                      if (j3 == j5 && j4 == j6) dshift->element(j3i, j6i) -= z * 2.0;
                      if (j3 == j5) dshift->element(j3i, j6i) += z * rdm1->element(j4, j6);
                      if (j4 == j6) dshift->element(j3i, j6i) += z * rdm1->element(j5, j3);
                      if (j4 == j5) dshift->element(j3i, j6i) -= z * 2.0 * rdm1->element(j3, j6);
                    }
                  }
                  dshift->element(j3i, j6i) += z * rdm1->element(j3, j6);
                }
              }
            }
          }
      ioffset += size_airj;
    }

    // r i s j
    {
      const size_t interm_size = denom_->shalf_hh()->ndim();
      for (size_t j5 = 0; j5 != nclo; ++j5)
        for (size_t j6 = 0; j6 != nclo; ++j6) {
          const size_t j5i = j5;
          const size_t j6i = j6;
          for (size_t j13 = 0; j13 != interm_size; ++j13) {
            const size_t jall = j5 + nclo * (j6 + nclo * j13) + ioffset;
            const double denom = - eig_[j5+ncore] - eig_[j6+ncore] + denom_->denom_hh(j13) - e0_;
            const double z = (*l)[jall] * (*t)[jall] * shift2 / (denom * denom);
            dshift->element(j5i, j5i) += z;
            dshift->element(j6i, j6i) += z;
            for (int j2 = 0; j2 != nact; ++j2) {
              const size_t j2i = j2 + nclo;
              for (int j3 = 0; j3 != nact; ++j3) {
                const size_t j3i = j3 + nclo;
                for (int j4 = 0; j4 != nact; ++j4) {
                  const size_t j4i = j4 + nclo;
                  dshift->element(j2i, j3i) -= 4.0 * z * denom_->shalf_hh()->element(j13, j4+j2*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j3*nact + istate*nact*nact);
                  dshift->element(j2i, j3i) += 2.0 * z * denom_->shalf_hh()->element(j13, j4+j2*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j4*nact + istate*nact*nact);
                  dshift->element(j4i, j3i) += 2.0 * z * denom_->shalf_hh()->element(j13, j4+j2*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j2+j3*nact + istate*nact*nact);
                  dshift->element(j4i, j3i) -= 4.0 * z * denom_->shalf_hh()->element(j13, j4+j2*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j2*nact + istate*nact*nact);
                  for (int j1 = 0; j1 != nact; ++j1) {
                    const size_t j1i = j1 + nclo;
                    dshift->element(j2i, j3i) -=       z * rdm1->element(j2, j1) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j3*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j2, j1) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j4*nact + istate*nact*nact);
                    dshift->element(j4i, j3i) -=       z * rdm1->element(j2, j1) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j2*nact + istate*nact*nact);
                    dshift->element(j4i, j3i) += 2.0 * z * rdm1->element(j2, j1) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j2+j3*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j2, j4) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j3*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j4, j3) * denom_->shalf_hh()->element(j13, j1+j2*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) -=       z * rdm1->element(j2, j1) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j4*nact + istate*nact*nact);
                    dshift->element(j1i, j3i) -=       z * rdm1->element(j2, j3) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j2+j4*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j4, j1) * denom_->shalf_hh()->element(j13, j2+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j4*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) -=       z * rdm1->element(j4, j1) * denom_->shalf_hh()->element(j13, j2+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j3*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) -=       z * rdm1->element(j4, j3) * denom_->shalf_hh()->element(j13, j2+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j4, j3) * denom_->shalf_hh()->element(j13, j2+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) -= 4.0 * z * rdm1->element(j2, j3) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact);
                    dshift->element(j2i, j3i) += 2.0 * z * rdm1->element(j2, j3) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact);
                    for (int j0 = 0; j0 != nact; ++j0) {
                      dshift->element(j2i, j3i) -=       z * rdm2->element(j2, j0, j4, j1) * denom_->shalf_hh()->element(j13, j0+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j3+j4*nact + istate*nact*nact);
                      dshift->element(j2i, j3i) -=       z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j0+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j3*nact + istate*nact*nact);
                      dshift->element(j2i, j3i) -=       z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j0+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact);
                      dshift->element(j2i, j3i) += 2.0 * z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j0+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact);
                      dshift->element(j1i, j3i) -=       z * rdm2->element(j2, j0, j4, j3) * denom_->shalf_hh()->element(j13, j0+j1*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j2+j4*nact + istate*nact*nact);
                      dshift->element(j2i, j3i) += 2.0 * z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j1+j0*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j4*nact + istate*nact*nact);
                      dshift->element(j2i, j3i) -=       z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j1+j0*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j4+j1*nact + istate*nact*nact);
                      dshift->element(j1i, j3i) -=       z * rdm2->element(j4, j0, j2, j3) * denom_->shalf_hh()->element(j13, j1+j0*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j2+j4*nact + istate*nact*nact);
                      for (int j7 = 0; j7 != nact; ++j7) {
                        dshift->element(j2i, j3i) -= z * rdm3->element(j0, j1, j4, j7, j2, j3) * denom_->shalf_hh()->element(j13, j0+j4*nact + istate*nact*nact) * denom_->shalf_hh()->element(j13, j1+j7*nact + istate*nact*nact);
                      }
                    }
                  }
                }
                dshift->element(j2i, j3i) += z * rdm1->element(j2, j3);
              }
            }
          }
        }
      ioffset += size_risj;
    }

    // a i r s & a r s i
    ioffset += size_airs;
    // a r s t
    {
      const size_t interm_size = denom_->shalf_xxh()->ndim();
      for (size_t j1 = 0; j1 != nvirt; ++j1) {
        size_t j1i = j1 + nocc - ncore;
        for (size_t j023 = 0; j023 != interm_size; ++j023) {
          const size_t jall = j1 + nvirt * j023 + ioffset;
          const double denom = eig_[j1+nocc] + denom_->denom_xxh(j023) - e0_;
          const double z = (*l)[jall] * (*t)[jall] * shift2 / (denom * denom);
          dshift->element(j1i, j1i) -= z;
          for (size_t j0 = 0; j0 != nact; ++j0) {
            const size_t j0i = j0 + nclo;
            for (size_t j2 = 0; j2 != nact; ++j2) {
              const size_t j2i = j2 + nclo;
              for (size_t j5 = 0; j5 != nact; ++j5) {
                const double VtuwO = denom_->shalf_xxh()->element(j023, j5+nact*(j0+nact*j2) + istate*nact*nact*nact);
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  const size_t j3i = j3 + nclo;
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    for (size_t j6 = 0; j6 != nact; ++j6) {
                      const double VxyzO = denom_->shalf_xxh()->element(j023, j6+nact*(j4+nact*j3) + istate*nact*nact*nact);
                      const double factor = VtuwO * VxyzO;
                      dshift->element(j2i, j3i) -= z * rdm2->element(j0, j4, j5, j6) * factor;
                      for (size_t j7 = 0; j7 != nact; ++j7) {
                        const size_t j7i = j7 + nclo;
                        dshift->element(j7i, j3i) -= z * rdm3->element(j7, j4, j0, j2, j5, j6) * factor;
                        dshift->element(j2i, j7i) -= z * rdm3->element(j3, j4, j5, j6, j0, j7) * factor;
                        for (size_t j8 = 0; j8 != nact; ++j8) {
                          const size_t j8i = j8 + nclo;
                          dshift->element(j7i, j8i) -= z * rdm4->element(j0, j2, j3, j4, j5, j6, j7, j8) * factor;
                          if (j3 == j2) {
                            dshift->element(j7i, j8i) -= z * rdm3->element(j0, j4, j5, j6, j7, j8) * factor;
                          }
                        }
                      }
                    }
                  }
                }
              }
              dshift->element(j0i, j2i) += z * rdm1->element(j0, j2);
            }
          }
        }
      }
      ioffset += size_arst;
    }

    // r i s t
  }

  return dshift;
}

#endif
