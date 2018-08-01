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


tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>,shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> CASPT2::CASPT2::feed_rdm(const int ist, const int jst) {
  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  shared_ptr<RDM<3>> rdm3;
  shared_ptr<RDM<4>> rdm4;

  const size_t nact = info_->nact();

  // collect den1ci
  {
    vector<IndexRange> o = rdm1all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    auto d1 = make_shared<RDM<1>>(nact);
    for (auto& i1 : o[1].range())
      for (auto& i0 : o[0].range()) {
        auto input = rdm1all_->at(jst, ist)->get_block(i0, i1);
        for (size_t io1 = 0; io1 != i1.size(); ++io1)
          copy_n(&input[0 + i0.size() * io1], i0.size(), d1->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1));
      }
    rdm1 = d1->copy();
  }

  // collect den2ci
  {
    vector<IndexRange>o = rdm2all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    auto d2 = make_shared<RDM<2>>(nact);
    for (auto& i3 : o[3].range())
      for (auto& i2 : o[2].range())
        for (auto& i1 : o[1].range())
          for (auto& i0 : o[0].range()) {
            auto input = rdm2all_->at(jst, ist)->get_block(i0, i1, i2, i3);
            for (size_t io3 = 0; io3 != i3.size(); ++io3)
              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * io3))], i0.size(),
                         d2->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2, io3 + i3.offset() - off3));
          }
    rdm2 = d2->copy();
  }

  // collect den3ci
  {
    vector<IndexRange>o = rdm3all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    auto d3 = make_shared<RDM<3>>(nact);
    for (auto& i5 : o[5].range())
      for (auto& i4 : o[4].range())
        for (auto& i3 : o[3].range())
          for (auto& i2 : o[2].range())
            for (auto& i1 : o[1].range())
              for (auto& i0 : o[0].range()) {
                auto input = rdm3all_->at(jst, ist)->get_block(i0, i1, i2, i3, i4, i5);
                for (size_t io5 = 0; io5 != i5.size(); ++io5)
                  for (size_t io4 = 0; io4 != i4.size(); ++io4)
                    for (size_t io3 = 0; io3 != i3.size(); ++io3)
                      for (size_t io2 = 0; io2 != i2.size(); ++io2)
                        for (size_t io1 = 0; io1 != i1.size(); ++io1)
                          copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * io5))))],
                                 i0.size(), d3->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                 io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5));
              }
    rdm3 = d3->copy();
  }

  // collect den4ci
  {
    vector<IndexRange>o = rdm4all_->at(jst, ist)->indexrange();
    const int off0 = o[0].front().offset();
    const int off1 = o[1].front().offset();
    const int off2 = o[2].front().offset();
    const int off3 = o[3].front().offset();
    const int off4 = o[4].front().offset();
    const int off5 = o[5].front().offset();
    const int off6 = o[6].front().offset();
    const int off7 = o[7].front().offset();
    auto d4 = make_shared<RDM<4>>(nact);
    for (auto& i7 : o[7].range())
      for (auto& i6 : o[6].range())
        for (auto& i5 : o[5].range())
          for (auto& i4 : o[4].range())
            for (auto& i3 : o[3].range())
              for (auto& i2 : o[2].range())
                for (auto& i1 : o[1].range())
                  for (auto& i0 : o[0].range()) {
                    auto input = rdm4all_->at(jst, ist)->get_block(i0, i1, i2, i3, i4, i5, i6, i7);
                    for (size_t io7 = 0; io7 != i7.size(); ++io7)
                      for (size_t io6 = 0; io6 != i6.size(); ++io6)
                        for (size_t io5 = 0; io5 != i5.size(); ++io5)
                          for (size_t io4 = 0; io4 != i4.size(); ++io4)
                            for (size_t io3 = 0; io3 != i3.size(); ++io3)
                              for (size_t io2 = 0; io2 != i2.size(); ++io2)
                                for (size_t io1 = 0; io1 != i1.size(); ++io1)
                                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * io5 + i5.size() * (io6 + i6.size() * io7)))))],
                                         i0.size(), d4->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                         io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5, io6 + i6.offset() - off6, io7 + i7.offset() - off7));
              }
    rdm4 = d4->copy();
  }

  return tie(rdm1, rdm2, rdm3, rdm4);
}


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
      {
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
      {
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
      {
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
      {
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
      {
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
      {
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


shared_ptr<Matrix> CASPT2::CASPT2::make_d2_imag(vector<shared_ptr<VectorB>> lambda, vector<shared_ptr<VectorB>> amplitude, vector<shared_ptr<VectorB>> residual) const {
  shared_ptr<Matrix> dshift = den2_->clone();
  const size_t nact = info_->nact();
  const size_t nclosed = info_->nclosed();
  const size_t nvirt = info_->nvirt();
  const size_t nocc = nact + nclosed;
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;
//  const size_t norb = nocc + nvirt;
  const size_t size_aibj = nvirt * nvirt * nclo * nclo;
  const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
  const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
  const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
  const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
  const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
  const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
  const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;
//  const double ovlfactor = 1.0;     //temporary

//  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;
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
    shared_ptr<const VectorB> r = residual[istate];
    size_t ioffset = 0;
    // a i b j
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
              const double Lambda = lcovar * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j0i, j0i) += Lambda;
              dshift->element(j1i, j1i) -= Lambda;
              dshift->element(j2i, j2i) += Lambda;
              dshift->element(j3i, j3i) -= Lambda;
            }
    }

    ioffset += size_aibj;
    ioffset += size_arbs;
    // a r b i
    {
      const size_t interm_size = denom_->shalf_x()->ndim();
      auto smallz = make_shared<Matrix>(interm_size, interm_size);
      auto largey = make_shared<Matrix>(interm_size, interm_size);
      auto largex = make_shared<Matrix>(interm_size, interm_size);
      auto dshift_part1 = make_shared<Matrix>(nact, nact);
      // multipliers
      for (size_t j3 = 0; j3 != nvirt; ++j3) {
        const size_t j3i = j3 + nocc - ncore;
        for (size_t j2 = 0; j2 != nclo; ++j2) {
          const size_t j2i = j2;
          for (size_t j1 = 0; j1 != nvirt; ++j1) {
            const size_t j1i = j1 + nocc - ncore;
            for (size_t j0o = 0; j0o != interm_size; ++j0o) {
              const size_t jall = j0o + interm_size * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
              const size_t jall2 = j0o + interm_size * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j0o) - e0all_[istate];
              const double Lambda = lcovar * (*t)[jall] * shift2 / (denom * denom);
              dshift->element(j1i, j1i) -= Lambda;
              dshift->element(j3i, j3i) -= Lambda;
              dshift->element(j2i, j2i) += Lambda;
              // should be done only once!!
              for (size_t j0 = 0; j0 != nact; ++j0) {
                const size_t j0i = j0 + nclo;
                for (size_t j6 = 0; j6 != nact; ++j6) {
                  const size_t j6i = j6 + nclo;
                  dshift->element(j0i, j6i) += Lambda * rdm1->element(j0, j6);
                  dshift_part1->element(j0, j6) += Lambda * rdm1->element(j0, j6);
                }
              }
              largey->element(j0o, j0o) -= Lambda * denom_->denom_x(j0o);
              for (size_t j1o = 0; j1o != interm_size; ++j1o) {
                const size_t kall = j1o + interm_size * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
                const double denom2 = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j1o) - e0all_[istate];
                largey->element(j0o, j1o) += lcovar * (*t)[kall] * shift2 * (1.0 / denom - 1.0 / denom2);
              }
#if 1
              for (size_t is = 0; is != nstates_; ++is) {
                for (size_t js = 0; js != nstates_; ++js) {
                  if (is != js) continue;
                  shared_ptr<const RDM<1>> rdm1tmp;
                  shared_ptr<const RDM<2>> rdm2tmp;
                  tie(rdm1tmp, rdm2tmp) = info_->rdm12(is, js);
                  for (size_t j0 = 0; j0 != nact; ++j0) {
                    for (size_t j6 = 0; j6 != nact; ++j6) {
                      const double VtO = denom_->shalf_x()->element(j0o, j0 + is*nact);
                      const double VuO = denom_->shalf_x()->element(j0o, j6 + js*nact);
                      const size_t j6i = j6 + nclo;
                      const size_t j0i = j0 + nclo;
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        for (size_t j5 = 0; j5 != nact; ++j5) {
                          dshift->element(j0i, j6i) -= Lambda * VtO * VuO * rdm2tmp->element(j0, j6, j4, j5);
                          dshift_part1->element(j0, j6) -= Lambda * VtO * VuO * rdm2tmp->element(j0, j6, j4, j5);
                        }
                      }
                    }
                  }
                }
              }
#endif
            }
          }
        }
      }
      fockact_->print("fockact = ");
      dshift_part1->print("dshift_part1 = ");
      dshift_part1->zero();
      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
        for (size_t j1o = 0; j1o != interm_size; ++j1o) {
          const double fdiff = denom_->denom_x(j1o) - denom_->denom_x(j0o);
          smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-8 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
        }
      }
      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
        for (size_t j1o = 0; j1o != interm_size; ++j1o) {
          largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_x(j1o))
                                    + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_x(j0o));
        }
      }
      smallz->print("z = ");
      largey->print("Y = ");
      largex->print("X = ");
      double xdiag = 0.0;
      for (size_t j0o = 0; j0o != interm_size; ++j0o)
        xdiag += largex->element(j0o, j0o);

#if 1
      auto f0_test = make_shared<Matrix>(interm_size, interm_size);
      for (size_t is = 0; is != nstates_; ++is) {
        for (size_t js = 0; js != nstates_; ++js) {
          if (is != js) continue;
          shared_ptr<const RDM<1>> rdm1tmp;
          shared_ptr<const RDM<2>> rdm2tmp;
          tie(rdm1tmp, rdm2tmp) = info_->rdm12(is, js);
          for (size_t j0 = 0; j0 != nact; ++j0) {
            const size_t j0i = j0 + nclo;
            for (size_t j6 = 0; j6 != nact; ++j6) {
              const size_t j6i = j6 + nclo;
              for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                for (size_t j1o = 0; j1o != interm_size; ++j1o) {
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      const double VtO = denom_->shalf_x()->element(j0o, j4 + is*nact);
                      const double VuO = denom_->shalf_x()->element(j1o, j5 + js*nact);
                      const double factor = VtO * VuO;
                      // How is fock? test.
                      f0_test->element(j0o, j1o) += VtO * VuO * rdm2tmp->element(j4, j5, j0, j6) * fockact_->element(j0, j6);
                      dshift->element(j0i, j6i) += smallz->element(j0o, j1o) * factor * rdm2tmp->element(j4, j5, j0, j6);
                      dshift_part1->element(j0, j6) += smallz->element(j0o, j1o) * factor * rdm2tmp->element(j4, j5, j0, j6);
                    }
                  }
                }
              }
            }
          }
        }
      }
      f0_test->print("f0 = ");
      for (size_t i = 0; i != interm_size; ++i) cout << setw(20) << setprecision(8) << denom_->denom_x(i) << endl;
      dshift_part1->print("dshift_part2 = ");
#endif
      double tes = 0.0;
      auto e1_test = make_shared<Matrix>(interm_size, interm_size);
      for (size_t is = 0; is != nstates_; ++is) {
        for (size_t js = 0; js != nstates_; ++js) {
          if (is != js) continue;
          shared_ptr<const RDM<1>> rdm1tmp;
          shared_ptr<const RDM<2>> rdm2tmp;
          tie(rdm1tmp, rdm2tmp) = info_->rdm12(is, js);
          auto e1 = make_shared<Matrix>(nact, nact);       // only for e^(1). we should implement RDM form of e other than e^(1).
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j6 = 0; j6 != nact; ++j6) {
              for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                for (size_t j1o = 0; j1o != interm_size; ++j1o) {
                  const double VtO = denom_->shalf_x()->element(j0o, j0 + is * nact);
                  const double VuO = denom_->shalf_x()->element(j1o, j6 + js * nact);
                  e1->element(j0, j6) -= largex->element(j0o, j1o) * VtO * VuO;
                  e1_test->element(j0o, j1o) += VtO * VuO * rdm1tmp->element(j0, j6);
                }
              }
              tes += e1->element(j0, j6) * rdm1tmp->element(j0, j6);
            }
          }
          e1->print("e1 = ");
        }
      }
      e1_test->print("e1_test = ");
      cout << " tes = " << tes << endl;
      cout << " xdiag = " << xdiag << endl;
    }
    ioffset += size_arbi;
    ioffset += size_airj;
    ioffset += size_risj;
    ioffset += size_airs;
    ioffset += size_arst;
    ioffset += size_rist;
  }

  return dshift;
}


double CASPT2::CASPT2::compute_energy_lt() {
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
  double E_lt = 0.0;

//  const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;
  const double shift2 = info_->shift() * info_->shift();
  for (int istate = 0; istate != nstates_; ++istate) { // state of T
    // temporary. will replace (RDM calculation redundant)
    // considering only SS-SR case
    shared_ptr<const VectorB> t = t2all_orthogonal()[istate];
    shared_ptr<const VectorB> l = lall_orthogonal()[istate];
    // a i b j
    size_t ioffset = 0;
    {
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j2 = 0; j2 != nclo; ++j2)
          for (size_t j1 = 0; j1 != nvirt; ++j1)
            for (size_t j0 = 0; j0 != nclo; ++j0) {
              const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
              const size_t jall2 = j0 + nclo * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
              const double lcovar = ((*l)[jall] * 8.0 - (*l)[jall2] * 4.0);
              const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j1+nocc] + eig_[j3+nocc];
              E_lt += lcovar * (*t)[jall] * shift2 / denom * 2.0;
            }
      ioffset += size_aibj;
    }

    // a r b s
    {
      const size_t interm_size = denom_->shalf_xx()->ndim();
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j1 = 0; j1 != nvirt; ++j1) {
            for (size_t j02 = 0; j02 != interm_size; ++j02) {
              const size_t jall = j02 + interm_size * (j1 + nvirt * j3) + ioffset;
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j02) - e0all_[istate];
              E_lt += (*l)[jall] * (*t)[jall] * shift2 / denom * 2.0;
            }
        }
    }
    ioffset += size_arbs;

    // a r b i
    {
      const size_t interm_size = denom_->shalf_x()->ndim();
      for (size_t j3 = 0; j3 != nvirt; ++j3)
        for (size_t j2 = 0; j2 != nclo; ++j2)
          for (size_t j1 = 0; j1 != nvirt; ++j1) {
            for (size_t j0o = 0; j0o != interm_size; ++j0o) {
              const size_t jall = j0o + interm_size * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
              const size_t jall2 = j0o + interm_size * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j0o) - e0all_[istate];
              E_lt += lcovar * (*t)[jall] * shift2 / denom * 2.0;
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
            for (size_t j3o = 0; j3o != interm_size; ++j3o) {
              const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3o)) + ioffset;
              const size_t jall2 = j2 + nclo * (j1 + nvirt * (j0 + nclo * j3o)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j3o) - e0all_[istate];
              E_lt += lcovar * (*t)[jall] * shift2 / denom * 2.0;
            }
          }
      ioffset += size_airj;
    }

    // r i s j
    {
      const size_t interm_size = denom_->shalf_hh()->ndim();
      for (size_t j5 = 0; j5 != nclo; ++j5)
        for (size_t j6 = 0; j6 != nclo; ++j6) {
          for (size_t j13 = 0; j13 != interm_size; ++j13) {
            const size_t jall = j5 + nclo * (j6 + nclo * j13) + ioffset;
            const double denom = - eig_[j5+ncore] - eig_[j6+ncore] + denom_->denom_hh(j13) - e0all_[istate];
            E_lt += (*l)[jall] * (*t)[jall] * shift2 / denom * 2.0;
          }
        }
      ioffset += size_risj;
    }

    // a i r s & a r s i
    {
      const size_t interm_size = denom_->shalf_xh()->ndim();
      for (size_t j6 = 0; j6 != nclo; ++j6) {
        for (size_t j7 = 0; j7 != nvirt; ++j7) {
          for (size_t jo = 0; jo != interm_size; ++jo) {
            const size_t jall = j6 + nclo * (j7 + nvirt * jo) + ioffset;
            const double denom = eig_[j7+nocc] - eig_[j6+ncore] + denom_->denom_xh(jo) - e0all_[istate];
            E_lt += (*l)[jall] * (*t)[jall] * shift2 / denom * 2.0;
          }
        }
      }
    }
    ioffset += size_airs;

    // a r s t
    {
      const size_t interm_size = denom_->shalf_xxh()->ndim();
      for (size_t j8 = 0; j8 != nvirt; ++j8) {
        for (size_t jo = 0; jo != interm_size; ++jo) {
          const size_t jall = j8 + nvirt * jo + ioffset;
          const double denom = eig_[j8+nocc] + denom_->denom_xxh(jo) - e0all_[istate];
          E_lt += (*l)[jall] * (*t)[jall] * shift2 / denom * 2.0;
        }
      }
      ioffset += size_arst;
    }

    // r i s t
    {
      const size_t interm_size = denom_->shalf_xhh()->ndim();
      for (size_t j8 = 0; j8 != nclo; ++j8) {
        for (size_t jo = 0; jo != interm_size; ++jo) {
          const size_t jall = j8 + nclo * jo + ioffset;
          const double denom = - eig_[j8+ncore] + denom_->denom_xhh(jo) - e0all_[istate];
          E_lt += (*l)[jall] * (*t)[jall] * shift2 / denom * 2.0;
        }
      }
    }
    ioffset += size_rist;
  }
 

  return E_lt;
}

#endif
