//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_shift.cc
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

#include <src/smith/caspt2/MSCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>,shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> MSCASPT2::MSCASPT2::feed_rdm(const int ist, const int jst) const {
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
                                  copy_n(&input[0 + i0.size() * (io1 + i1.size() * (io2 + i2.size() * (io3 + i3.size() * (io4 + i4.size() * (io5 + i5.size() * (io6 + i6.size() * io7))))))],
                                         i0.size(), d4->element_ptr(i0.offset() - off0, io1 + i1.offset() - off1, io2 + i2.offset() - off2,
                                         io3 + i3.offset() - off3, io4 + i4.offset() - off4, io5 + i5.offset() - off5, io6 + i6.offset() - off6, io7 + i7.offset() - off7));
              }
    rdm4 = d4->copy();
  }

  return tie(rdm1, rdm2, rdm3, rdm4);
}


tuple<shared_ptr<Matrix>,shared_ptr<Vec<double>>,shared_ptr<VecRDM<1>>,shared_ptr<VecRDM<2>>,shared_ptr<VecRDM<3>>,shared_ptr<VecRDM<3>>,vector<double>>
  MSCASPT2::MSCASPT2::make_d2_imag(vector<shared_ptr<VectorB>> lambda, vector<shared_ptr<VectorB>> amplitude) const {
  const size_t nact = info_->nact();
  const int nstates = nact ? info_->ciwfn()->nstates() : 1;
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

  Timer timer(1);

  shared_ptr<Matrix> dshift = den2_->clone();
  auto e0 = make_shared<Vec<double>>();
  auto e1 = make_shared<VecRDM<1>>();
  auto e2 = make_shared<VecRDM<2>>();
  auto e3 = make_shared<VecRDM<3>>();
  auto e4 = make_shared<VecRDM<3>>();
  vector<double> nimag;
  nimag.resize(nstates);

  for (int i = 0; i != nstates; ++i)
    nimag[i] = 0.0;

  if (nact) {
    for (size_t is = 0; is != nstates; ++is)
      for (size_t js = 0; js != nstates; ++js)
        if (!info_->sssr() || is == js) {
          auto e0temp = make_shared<double>(0.0);
          auto e1temp = make_shared<RDM<1>>(nact);
          auto e2temp = make_shared<RDM<2>>(nact);
          auto e3temp = make_shared<RDM<3>>(nact);
          auto e4temp = make_shared<RDM<3>>(nact);

          e0->emplace(is, js, e0temp);
          e1->emplace(is, js, e1temp);
          e2->emplace(is, js, e2temp);
          e3->emplace(is, js, e3temp);
          e4->emplace(is, js, e4temp);
        }
  }

  const double shift2 = info_->shift() * info_->shift();
  size_t ioffset = 0;

  if (size_arbs) {
    const size_t interm_size = denom_->shalf_xx()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j3 = 0; j3 != nvirt; ++j3) {
        const size_t j3i = j3 + nocc - ncore;
        for (size_t j1 = 0; j1 != nvirt; ++j1) {
          const size_t j1i = j1 + nocc - ncore;
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            const size_t jall = j0o + interm_size * (j1 + nvirt * j3) + ioffset;
            const double denom = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j0o) - e0all_[istate];
            const double Lambda = shift2 * (*l)[jall] * (*t)[jall];
            dshift->element(j1i, j1i) -= Lambda;
            dshift->element(j3i, j3i) -= Lambda;
            smallz->element(j0o, j0o) -= Lambda;
            nimag[istate] -= Lambda;
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              const size_t kall = j1o + interm_size * (j1 + nvirt * j3) + ioffset;
              const double denomk = eig_[j3+nocc] + eig_[j1+nocc] + denom_->denom_xx(j1o) - e0all_[istate];
              const double lt = (*l)[jall] * (*t)[kall] * shift2 * denom;
              const double tl = (*t)[jall] * (*l)[kall] * shift2 * denomk;
              largey->element(j0o, j1o) += (*l)[jall] * (*t)[kall] * shift2 * (denom - denomk);
              largeq->element(j0o, j1o) += (lt + tl);
            }
          }
        }
      }
    }

    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_xx(j1o) - denom_->denom_xx(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_xx(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_xx(j0o));
      }
    }
    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_xx()->get_submatrix(0, is * nact * nact, interm_size, nact * nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<2>>(nact);
          // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j2 = 0; j2 != nact; ++j2) {
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                    const double VrsO = denom_->shalf_xx()->element(j0o, j2 + j3 * nact + js * nact * nact);
                    Rmat->element(j0, j1, j2, j3) += Qmat->element(j0o, j0 + j1 * nact) * VrsO;
                    e2->at(is, js)->element(j0, j2, j1, j3) += Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                  }
                }
              }
            }
          }
          // (3) form d^{(2)}_{rs} and e3
          for (size_t j0 = 0; j0 != nact; ++j0)
            for (size_t j1 = 0; j1 != nact; ++j1)
              for (size_t j2 = 0; j2 != nact; ++j2)
                for (size_t j3 = 0; j3 != nact; ++j3)
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    const size_t j4i = j4 + nclo;
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      const size_t j5i = j5 + nclo;
                      dshift->element(j4i, j5i) += Rmat->element(j0, j1, j2, j3) * rdm3->element(j0, j2, j1, j3, j4, j5);
                      e3->at(is,js)->element(j0, j2, j1, j3, j4, j5) += 2.0 * Rmat->element(j0, j1, j2, j3) * fockact_->element(j4, j5);
                    }
                  }
        }
      }
    }
    ioffset += size_arbs;
    timer.tick_print("dshift arbs");
  }

    // a r b i
#if 1
  if (size_arbi) {
    const size_t interm_size = denom_->shalf_x()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

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
              const double tcovar = ((*t)[jall] * 2.0 - (*t)[jall2]);
              const double denom = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j0o) - e0all_[istate];
              const double Lambda = shift2 * lcovar * (*t)[jall];
              dshift->element(j1i, j1i) -= Lambda;
              dshift->element(j3i, j3i) -= Lambda;
              dshift->element(j2i, j2i) += Lambda;
              smallz->element(j0o, j0o) -= Lambda;
              nimag[istate] -= Lambda;
              for (size_t j1o = 0; j1o != interm_size; ++j1o) {
                const size_t kall = j1o + interm_size * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
                const double denomk = eig_[j3+nocc] + eig_[j1+nocc] - eig_[j2+ncore] + denom_->denom_x(j1o) - e0all_[istate];
                const double lt = lcovar * (*t)[kall] * shift2 * denom;
                const double tl = tcovar * (*l)[kall] * shift2 * denomk;
                largey->element(j0o, j1o) += lcovar * (*t)[kall] * shift2 * (denom - denomk);
                largeq->element(j0o, j1o) += (lt + tl);
              }
            }
          }
        }
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_x(j1o) - denom_->denom_x(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_x(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_x(j0o));
      }
    }

    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_x()->get_submatrix(0, is * nact, interm_size, nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<1>>(nact);
          // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                const double VrO = denom_->shalf_x()->element(j0o, j1 + js * nact);
                Rmat->element(j0, j1) += Qmat->element(j0o, j0) * VrO;
                e1->at(is, js)->element(j0, j1) += Pmat->element(j0o, j0) * VrO;
              }
            }
          }
          // (3) form d^{(2)}_{rs} and e3
          for (size_t j0 = 0; j0 != nact; ++j0)
            for (size_t j1 = 0; j1 != nact; ++j1)
              for (size_t j2 = 0; j2 != nact; ++j2) {
                const size_t j2i = j2 + nclo;
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  const size_t j3i = j3 + nclo;
                  dshift->element(j2i, j3i) += Rmat->element(j0, j1) * rdm2->element(j0, j1, j2, j3);
                  e2->at(is, js)->element(j0, j1, j2, j3) += 2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3);
                }
              }
        }
      }
    }
    ioffset += size_arbi;
    timer.tick_print("dshift arbi");
  }
#endif

    // a i r j
  if (size_airj) {
    const size_t interm_size = denom_->shalf_h()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j2 = 0; j2 != nclo; ++j2) {
        const size_t j2i = j2;
        for (size_t j1 = 0; j1 != nvirt; ++j1) {
          const size_t j1i = j1 + nocc - ncore;
          for (size_t j0 = 0; j0 != nclo; ++j0) {
            const size_t j0i = j0;
            for (size_t j0o = 0; j0o != interm_size; ++j0o) {
              const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j0o)) + ioffset;
              const size_t jall2 = j2 + nclo * (j1 + nvirt * (j0 + nclo * j0o)) + ioffset;
              const double lcovar = ((*l)[jall] * 2.0 - (*l)[jall2]);
              const double tcovar = ((*t)[jall] * 2.0 - (*t)[jall2]);
              const double denom = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j0o) - e0all_[istate];
              const double Lambda = shift2 * lcovar * (*t)[jall];
              dshift->element(j1i, j1i) -= Lambda;
              dshift->element(j0i, j0i) += Lambda;
              dshift->element(j2i, j2i) += Lambda;
              smallz->element(j0o, j0o) -= Lambda;
              nimag[istate] -= Lambda;
              for (size_t j1o = 0; j1o != interm_size; ++j1o) {
                const size_t kall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j1o)) + ioffset;
                const double denomk = eig_[j1+nocc] - eig_[j0+ncore] - eig_[j2+ncore] + denom_->denom_h(j1o) - e0all_[istate];
                const double lt = lcovar * (*t)[kall] * shift2 * denom;
                const double tl = tcovar * (*l)[kall] * shift2 * denomk;
                largey->element(j0o, j1o) += lcovar * (*t)[kall] * shift2 * (denom - denomk);
                largeq->element(j0o, j1o) += (lt + tl);
              }
            }
          }
        }
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_h(j1o) - denom_->denom_h(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_h(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_h(j0o));
      }
    }

    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_h()->get_submatrix(0, is * nact, interm_size, nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<1>>(nact);
          // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                const double VrO = denom_->shalf_h()->element(j0o, j1 + js * nact);
                Rmat->element(j0, j1) += Qmat->element(j0o, j0) * VrO;
                e1->at(is, js)->element(j0, j1) += -Pmat->element(j0o, j0) * VrO;
                (*e0->at(is, js)) += Pmat->element(j0o, j0) * VrO * 2.0;
              }
            }
          }
          // (3) form d^{(2)}_{rs} and e3
          for (size_t j0 = 0; j0 != nact; ++j0)
            for (size_t j1 = 0; j1 != nact; ++j1)
              for (size_t j2 = 0; j2 != nact; ++j2) {
                const size_t j2i = j2 + nclo;
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  const size_t j3i = j3 + nclo;
                  dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm2->element(j0, j1, j2, j3);
                  if (j1 == j2 && j0 == j3 && is == js) dshift->element(j2i, j3i) +=  2.0 * Rmat->element(j0, j1);
                  if (j1 == j2)                         dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm1->element(j0, j3);
                  if (j0 == j3)                         dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm1->element(j1, j2);
                  if (j0 == j1 && is == js)             dshift->element(j2i, j3i) +=  2.0 * Rmat->element(j0, j1) * rdm1->element(j2, j3);
                  e2->at(is, js)->element(j0, j1, j2, j3) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                  if (j1 == j2 && j0 == j3 && is == js) (*e0->at(is, js)) += 2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                  if (j1 == j2)                         e1->at(is, js)->element(j0, j3) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                  if (j0 == j3)                         e1->at(is, js)->element(j1, j2) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                  if (j0 == j1 && is == js)             e1->at(is, js)->element(j2, j3) +=  2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                }
              }
        }
      }
    }
    ioffset += size_airj;
    timer.tick_print("dshift airj");
  }

    // r i s j
  if (size_risj) {
    const size_t interm_size = denom_->shalf_hh()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j0 = 0; j0 != nclo; ++j0) {
        const size_t j0i = j0;
        for (size_t j1 = 0; j1 != nclo; ++j1) {
          const size_t j1i = j1;
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            const size_t jall = j0 + nclo * (j1 + nclo * j0o) + ioffset;
            const double denom = - eig_[j0+ncore] - eig_[j1+ncore] + denom_->denom_hh(j0o) - e0all_[istate];
            const double Lambda = shift2 * (*l)[jall] * (*t)[jall];
            dshift->element(j0i, j0i) += Lambda;
            dshift->element(j1i, j1i) += Lambda;
            smallz->element(j0o, j0o) -= Lambda;
            nimag[istate] -= Lambda;
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              const size_t kall = j0 + nclo * (j1 + nclo * j1o) + ioffset;
              const double denomk = - eig_[j0+ncore] - eig_[j1+ncore] + denom_->denom_hh(j1o) - e0all_[istate];
              const double lt = (*l)[jall] * (*t)[kall] * shift2 * denom;
              const double tl = (*t)[jall] * (*l)[kall] * shift2 * denomk;
              largey->element(j0o, j1o) += (*l)[jall] * (*t)[kall] * shift2 * (denom - denomk);
              largeq->element(j0o, j1o) += (lt + tl);
            }
          }
        }
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_hh(j1o) - denom_->denom_hh(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_hh(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_hh(j0o));
      }
    }

    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_hh()->get_submatrix(0, is * nact * nact, interm_size, nact * nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<2>>(nact);
          // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j2 = 0; j2 != nact; ++j2) {
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                    const double VrsO = denom_->shalf_hh()->element(j0o, j2 + j3 * nact + js * nact * nact);
                    Rmat->element(j0, j2, j1, j3) += Qmat->element(j0o, j0 + j1 * nact) * VrsO;
                    e2->at(is, js)->element(j0, j2, j1, j3) += Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j2 == j1)                         e1->at(is, js)->element(j0, j3) +=  1.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j1 == j3)                         e1->at(is, js)->element(j0, j2) += -2.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j0 == j2)                         e1->at(is, js)->element(j1, j3) += -2.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j0 == j3)                         e1->at(is, js)->element(j1, j2) +=  1.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j0 == j2 && j1 == j3 && is == js) *(e0->at(is, js)) +=  4.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                    if (j0 == j3 && j1 == j2 && is == js) *(e0->at(is, js)) += -2.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                  }
                }
              }
            }
          }
          // (3) form d^{(2)}_{rs} and e3
          for (size_t j4 = 0; j4 != nact; ++j4)
            for (size_t j1 = 0; j1 != nact; ++j1)
              for (size_t j5 = 0; j5 != nact; ++j5)
                for (size_t j0 = 0; j0 != nact; ++j0)
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    const size_t j2i = j2 + nclo;
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      const size_t j3i = j3 + nclo;
                      const double factor = Rmat->element(j0, j5, j1, j4);
                      dshift->element(j2i, j3i) += factor * rdm3->element(j0, j5, j1, j4, j2, j3);
                      if (j2 == j5)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j1, j4, j0, j3);
                      if (j2 == j4)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j0, j5, j1, j3);
                      if (j1 == j5)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j0, j4, j2, j3);
                      if (j1 == j5 && j2 == j4)                         dshift->element(j2i, j3i) += factor * rdm1->element(j0, j3);
                      if (j1 == j4)                                     dshift->element(j2i, j3i) += factor * -2.0 * rdm2->element(j0, j5, j2, j3);
                      if (j1 == j4 && j2 == j5)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j0, j3);
                      if (j1 == j3)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j0, j5, j2, j4);
                      if (j1 == j3 && j2 == j5)                         dshift->element(j2i, j3i) += factor * rdm1->element(j0, j4);
                      if (j1 == j3 && j2 == j4)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j0, j5);
                      if (j0 == j5)                                     dshift->element(j2i, j3i) += factor * -2.0 * rdm2->element(j1, j4, j2, j3);
                      if (j0 == j5 && j2 == j4)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j1, j3);
                      if (j1 == j3 && j0 == j5)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j2, j4);
                      if (j1 == j3 && j0 == j5 && j2 == j4 && is == js) dshift->element(j2i, j3i) += factor *  4.0;
                      if (j0 == j4)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j1, j5, j2, j3);
                      if (j0 == j4 && j2 == j5)                         dshift->element(j2i, j3i) += factor * rdm1->element(j1, j3);
                      if (j1 == j3 && j0 == j4)                         dshift->element(j2i, j3i) += factor * rdm1->element(j2, j5);
                      if (j1 == j3 && j0 == j4 && j2 == j5 && is == js) dshift->element(j2i, j3i) += factor * -2.0;
                      if (j0 == j3)                                     dshift->element(j2i, j3i) += factor * rdm2->element(j2, j5, j1, j4);
                      if (j0 == j3 && j2 == j5)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j1, j4);
                      if (j0 == j3 && j2 == j4)                         dshift->element(j2i, j3i) += factor * rdm1->element(j1, j5);
                      if (j1 == j5 && j0 == j3)                         dshift->element(j2i, j3i) += factor * rdm1->element(j2, j4);
                      if (j0 == j3 && j1 == j5 && j2 == j4 && is == js) dshift->element(j2i, j3i) += factor * -2.0;
                      if (j0 == j3 && j1 == j4)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j2, j5);
                      if (j1 == j4 && j2 == j5 && j0 == j3 && is == js) dshift->element(j2i, j3i) += factor *  4.0;
                      if (j0 == j5 && j1 == j4)                         dshift->element(j2i, j3i) += factor *  4.0 * rdm1->element(j2, j3);
                      if (j0 == j4 && j1 == j5)                         dshift->element(j2i, j3i) += factor * -2.0 * rdm1->element(j2, j3);
                      e3->at(is,js)->element(j0, j5, j1, j4, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j2 == j5)                                     e2->at(is,js)->element(j1, j4, j0, j3) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j2 == j4)                                     e2->at(is,js)->element(j0, j5, j1, j3) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j5)                                     e2->at(is,js)->element(j0, j4, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j5 && j2 == j4)                         e1->at(is,js)->element(j0, j3)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j4)                                     e2->at(is,js)->element(j0, j5, j2, j3) += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j4 && j2 == j5)                         e1->at(is,js)->element(j0, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3)                                     e2->at(is,js)->element(j0, j5, j2, j4) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j2 == j5)                         e1->at(is,js)->element(j0, j4)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j2 == j4)                         e1->at(is,js)->element(j0, j5)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j5)                                     e2->at(is,js)->element(j1, j4, j2, j3) += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j5 && j2 == j4)                         e1->at(is,js)->element(j1, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j0 == j5)                         e1->at(is,js)->element(j2, j4)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j0 == j5 && j2 == j4 && is == js) (*e0->at(is,js))                       += 2.0 *  4.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j4)                                     e2->at(is,js)->element(j1, j5, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j4 && j2 == j5)                         e1->at(is,js)->element(j1, j3)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j0 == j4)                         e1->at(is,js)->element(j2, j5)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j3 && j0 == j4 && j2 == j5 && is == js) (*e0->at(is,js))                       += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j3)                                     e2->at(is,js)->element(j2, j5, j1, j4) += 2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j3 && j2 == j5)                         e1->at(is,js)->element(j1, j4)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j3 && j2 == j4)                         e1->at(is,js)->element(j1, j5)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j5 && j0 == j3)                         e1->at(is,js)->element(j2, j4)         += 2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j3 && j1 == j5 && j2 == j4 && is == js) (*e0->at(is,js))                       += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j3 && j1 == j4)                         e1->at(is,js)->element(j2, j5)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                      if (j1 == j4 && j2 == j5 && j0 == j3 && is == js) (*e0->at(is,js))                       += 2.0 *  4.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j5 && j1 == j4)                         e1->at(is,js)->element(j2, j3)         += 2.0 *  4.0 * factor * fockact_->element(j2, j3);
                      if (j0 == j4 && j1 == j5)                         e1->at(is,js)->element(j2, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                    }
                  }
        }
      }
    }
    ioffset += size_risj;
    timer.tick_print("dshift risj");
  }

    // a i r s & a r s i
  if (size_airs) {
    const size_t interm_size = denom_->shalf_xh()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j0 = 0; j0 != nclo; ++j0) {
        const size_t j0i = j0;
        for (size_t j1 = 0; j1 != nvirt; ++j1) {
          const size_t j1i = j1 + nocc - ncore;
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            const size_t jall = j0 + nclo * (j1 + nvirt * j0o) + ioffset;
            const double denom = eig_[j1+nocc] - eig_[j0+ncore] + denom_->denom_xh(j0o) - e0all_[istate];
            const double Lambda = shift2 * (*l)[jall] * (*t)[jall];
            dshift->element(j0i, j0i) += Lambda;
            dshift->element(j1i, j1i) -= Lambda;
            smallz->element(j0o, j0o) -= Lambda;
            nimag[istate] -= Lambda;
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              const size_t kall = j0 + nclo * (j1 + nvirt * j1o) + ioffset;
              const double denomk = eig_[j1+nocc] - eig_[j0+ncore] + denom_->denom_xh(j1o) - e0all_[istate];
              const double lt = (*l)[jall] * (*t)[kall] * shift2 * denom;
              const double tl = (*t)[jall] * (*l)[kall] * shift2 * denomk;
              largey->element(j0o, j1o) += (*l)[jall] * (*t)[kall] * shift2 * (denom - denomk);
              largeq->element(j0o, j1o) += (lt + tl);
            }
          }
        }
      }
    }

    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_xh(j1o) - denom_->denom_xh(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_xh(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_xh(j0o));
      }
    }

    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> VmatA = denom_->shalf_xh()->get_submatrix(0, (2 * is + 0) * nact * nact, interm_size, nact * nact);
          shared_ptr<Matrix> VmatB = denom_->shalf_xh()->get_submatrix(0, (2 * is + 1) * nact * nact, interm_size, nact * nact);
          auto QmatA = make_shared<Matrix>(*smallz % (2.0 * (*VmatA) - (*VmatB)));
          auto QmatB = make_shared<Matrix>(*smallz % (-1.0 * (*VmatA)));
          auto QmatC = make_shared<Matrix>(*smallz % (*VmatB));
          auto PmatA = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (2.0 * (*VmatA) - (*VmatB)));
          auto PmatB = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (-1.0 * (*VmatA)));
          auto PmatC = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*VmatB));

          auto RmatA = make_shared<RDM<2>>(nact);
          auto RmatC = make_shared<RDM<2>>(nact);
          // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j2 = 0; j2 != nact; ++j2) {
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                    const double VrsO = denom_->shalf_xh()->element(j0o, j2 + j3 * nact + (2 * js + 0) * nact * nact);
                    const double VrsS = denom_->shalf_xh()->element(j0o, j2 + j3 * nact + (2 * js + 1) * nact * nact);
                    RmatA->element(j0, j1, j3, j2) += QmatA->element(j0o, j0 + j1 * nact) * VrsO + QmatB->element(j0o, j0 + j1 * nact) * VrsS;
                    RmatC->element(j0, j1, j3, j2) += QmatC->element(j0o, j0 + j1 * nact) * VrsS;
                    e2->at(is, js)->element(j0, j1, j3, j2)       += (PmatA->element(j0o, j0 + j1 * nact) * VrsO + PmatB->element(j0o, j0 + j1 * nact) * VrsS);
                    if (j1 == j3) e1->at(is, js)->element(j0, j2) += (PmatA->element(j0o, j0 + j1 * nact) * VrsO + PmatB->element(j0o, j0 + j1 * nact) * VrsS);
                    e2->at(is, js)->element(j1, j3, j2, j0)       += -1.0 * PmatC->element(j0o, j0 + j1 * nact) * VrsS;
                    if (j1 == j3) e1->at(is, js)->element(j2, j0) +=  2.0 * PmatC->element(j0o, j0 + j1 * nact) * VrsS;
                  }
                }
              }
            }
          }
          // (3) form d^{(2)}_{rs} and e3
          for (size_t j4 = 0; j4 != nact; ++j4)
            for (size_t j1 = 0; j1 != nact; ++j1)
              for (size_t j5 = 0; j5 != nact; ++j5)
                for (size_t j0 = 0; j0 != nact; ++j0)
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    const size_t j2i = j2 + nclo;
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      const size_t j3i = j3 + nclo;
                      dshift->element(j2i, j3i) += RmatA->element(j0, j1, j4, j5) * rdm3->element(j0, j1, j4, j5, j2, j3);
                      if (j3 == j4)             dshift->element(j2i, j3i) += RmatA->element(j0, j1, j4, j5) * rdm2->element(j0, j1, j2, j5);
                      if (j1 == j2)             dshift->element(j2i, j3i) += RmatA->element(j0, j1, j4, j5) * rdm2->element(j0, j3, j4, j5);
                      if (j1 == j2 && j3 == j4) dshift->element(j2i, j3i) += RmatA->element(j0, j1, j4, j5) * rdm1->element(j0, j5);
                      if (j1 == j4)             dshift->element(j2i, j3i) += RmatA->element(j0, j1, j4, j5) * rdm2->element(j2, j3, j0, j5);
                      dshift->element(j2i, j3i) -= RmatC->element(j0, j1, j4, j5) * rdm3->element(j1, j4, j5, j0, j2, j3);
                      if (j3 == j1)             dshift->element(j2i, j3i) += -1.0 * RmatC->element(j0, j1, j4, j5) * rdm2->element(j2, j4, j5, j0);
                      if (j4 == j2)             dshift->element(j2i, j3i) += -1.0 * RmatC->element(j0, j1, j4, j5) * rdm2->element(j1, j3, j5, j0);
                      if (j3 == j1 && j4 == j2) dshift->element(j2i, j3i) +=  2.0 * RmatC->element(j0, j1, j4, j5) * rdm1->element(j5, j0);
                      if (j4 == j1)             dshift->element(j2i, j3i) +=  2.0 * RmatC->element(j0, j1, j4, j5) * rdm2->element(j2, j3, j5, j0);
                      e3->at(is,js)->element(j0, j1, j4, j5, j2, j3) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j3 == j4)             e2->at(is,js)->element(j0, j1, j2, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j1 == j2)             e2->at(is,js)->element(j0, j3, j4, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j1 == j2 && j3 == j4) e1->at(is,js)->element(j0, j5)         += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j1 == j4)             e2->at(is,js)->element(j2, j3, j0, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      e3->at(is,js)->element(j1, j4, j5, j0, j2, j3) -= 2.0 * RmatC->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j3 == j1)             e2->at(is,js)->element(j2, j4, j5, j0) += 2.0 * -1.0 * RmatC->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j4 == j2)             e2->at(is,js)->element(j1, j3, j5, j0) += 2.0 * -1.0 * RmatC->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j3 == j1 && j4 == j2) e1->at(is,js)->element(j5, j0)         += 2.0 *  2.0 * RmatC->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                      if (j4 == j1)             e2->at(is,js)->element(j2, j3, j5, j0) += 2.0 *  2.0 * RmatC->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                    }
                  }
        }
      }
    }
    ioffset += size_airs;
    timer.tick_print("dshift airs");
  }

    // a r s t
  if (size_arst) {
    const size_t interm_size = denom_->shalf_xxh()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j0 = 0; j0 != nvirt; ++j0) {
        const size_t j0i = j0 + nocc - ncore;
        for (size_t j0o = 0; j0o != interm_size; ++j0o) {
          const size_t jall = j0 + nvirt * j0o + ioffset;
          const double denom = eig_[j0+nocc] + denom_->denom_xxh(j0o) - e0all_[istate];
          const double Lambda = shift2 * (*l)[jall] * (*t)[jall];
          dshift->element(j0i, j0i) -= Lambda;
          smallz->element(j0o, j0o) -= Lambda;
          nimag[istate] -= Lambda;
          for (size_t j1o = 0; j1o != interm_size; ++j1o) {
            const size_t kall = j0 + nvirt * j1o + ioffset;
            const double denomk = eig_[j0+nocc] + denom_->denom_xxh(j1o) - e0all_[istate];
            const double lt = (*l)[jall] * (*t)[kall] * shift2 * denom;
            const double tl = (*t)[jall] * (*l)[kall] * shift2 * denomk;
            largey->element(j0o, j1o) += (*l)[jall] * (*t)[kall] * shift2 * (denom - denomk);
            largeq->element(j0o, j1o) += (lt + tl);
          }
        }
      }
    }

    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_xxh(j1o) - denom_->denom_xxh(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_xxh(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_xxh(j0o));
      }
    }
    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_xxh()->get_submatrix(0, is * nact * nact * nact, interm_size, nact * nact * nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<3>>(nact);
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j2 = 0; j2 != nact; ++j2) {
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                        const double VrstO = denom_->shalf_xxh()->element(j0o, j3 + nact * (j4 + nact * j5) + js * nact * nact * nact);
                        Rmat->element(j1, j2, j5, j4, j0, j3) += Qmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        e3->at(is, js)->element(j0, j1, j5, j4, j2, j3) += Pmat->element(j0o, j2 + nact * (j0 + nact * j1)) * VrstO;
                        if (j1 == j5) e2->at(is, js)->element(j0, j4, j2, j3) += Pmat->element(j0o, j2 + nact * (j0 + nact * j1)) * VrstO;
                      }
                    }
                  }
                }
              }
            }
          }
          for (size_t j7 = 0; j7 != nact; ++j7)
            for (size_t j0 = 0; j0 != nact; ++j0)
              for (size_t j6 = 0; j6 != nact; ++j6)
                for (size_t j5 = 0; j5 != nact; ++j5)
                  for (size_t j2 = 0; j2 != nact; ++j2)
                    for (size_t j1 = 0; j1 != nact; ++j1) {
                      e4->at(is,js)->element(j1, j2, j5, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7);
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        size_t j4i = j4 + nclo;
                        for (size_t j3 = 0; j3 != nact; ++j3) {
                          size_t j3i = j3 + nclo;
                          dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm4->element(j1, j2, j5, j6, j0, j7, j3, j4);
                          if (j4 == j5)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j1, j2, j3, j6, j0, j7);
                          if (j2 == j3)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j1, j4, j5, j6, j0, j7);
                          if (j2 == j3 && j4 == j5) dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm2->element(j1, j6, j0, j7);
                          if (j2 == j5)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j3, j4, j1, j6, j0, j7);
                          if (j4 == j5)             e3->at(is,js)->element(j1, j2, j3, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                          if (j2 == j3)             e3->at(is,js)->element(j1, j4, j5, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                          if (j2 == j3 && j4 == j5) e2->at(is,js)->element(j1, j6, j0, j7)         += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                          if (j2 == j5)             e3->at(is,js)->element(j3, j4, j1, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                        }
                      }
                    }
        }
      }
    }
    ioffset += size_arst;
    timer.tick_print("dshift arst");
  }

    // r i s t
  if (size_rist) {
    const size_t interm_size = denom_->shalf_xhh()->ndim();
    auto smallz = make_shared<Matrix>(interm_size, interm_size);
    auto largey = make_shared<Matrix>(interm_size, interm_size);
    auto largex = make_shared<Matrix>(interm_size, interm_size);
    auto largeq = make_shared<Matrix>(interm_size, interm_size);
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];

      for (size_t j0 = 0; j0 != nclo; ++j0) {
        const size_t j0i = j0;
        for (size_t j0o = 0; j0o != interm_size; ++j0o) {
          const size_t jall = j0 + nclo * j0o + ioffset;
          const double denom = - eig_[j0+ncore] + denom_->denom_xhh(j0o) - e0all_[istate];
          const double Lambda = shift2 * (*l)[jall] * (*t)[jall];
          dshift->element(j0i, j0i) += Lambda;
          smallz->element(j0o, j0o) -= Lambda;
          nimag[istate] -= Lambda;
          for (size_t j1o = 0; j1o != interm_size; ++j1o) {
            const size_t kall = j0 + nclo * j1o + ioffset;
            const double denomk = - eig_[j0+ncore] + denom_->denom_xhh(j1o) - e0all_[istate];
            const double lt = (*l)[jall] * (*t)[kall] * shift2 * denom;
            const double tl = (*t)[jall] * (*l)[kall] * shift2 * denomk;
            largey->element(j0o, j1o) += (*l)[jall] * (*t)[kall] * shift2 * (denom - denomk);
            largeq->element(j0o, j1o) += (lt + tl);
          }
        }
      }
    }

    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        if (j1o == j0o) continue;
        const double fdiff = denom_->denom_xhh(j1o) - denom_->denom_xhh(j0o);
        smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
      }
    }
    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
      for (size_t j1o = 0; j1o != interm_size; ++j1o) {
        largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * denom_->denom_xhh(j1o))
                                  + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * denom_->denom_xhh(j0o));
      }
    }
    for (size_t is = 0; is != nstates; ++is) {
      for (size_t js = 0; js != nstates; ++js) {
        if (info_->sssr() && is != js) continue;
        shared_ptr<RDM<1>> rdm1;
        shared_ptr<RDM<2>> rdm2;
        shared_ptr<RDM<3>> rdm3;
        shared_ptr<RDM<4>> rdm4;
        tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);

        {
          shared_ptr<Matrix> Vmat = denom_->shalf_xhh()->get_submatrix(0, is * nact * nact * nact, interm_size, nact * nact * nact);
          auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
          auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));

          auto Rmat = make_shared<RDM<3>>(nact);
          for (size_t j0 = 0; j0 != nact; ++j0) {
            for (size_t j1 = 0; j1 != nact; ++j1) {
              for (size_t j2 = 0; j2 != nact; ++j2) {
                for (size_t j3 = 0; j3 != nact; ++j3) {
                  for (size_t j4 = 0; j4 != nact; ++j4) {
                    for (size_t j5 = 0; j5 != nact; ++j5) {
                      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                        const double VrstO = denom_->shalf_xhh()->element(j0o, j3 + nact * (j4 + nact * j5) + js * nact * nact * nact);
                        Rmat->element(j0, j1, j5, j2, j4, j3) += Qmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        e3->at(is,js)->element(j0, j1, j5, j2, j4, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j2 == j4)             e2->at(is, js)->element(j0, j1, j5, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j2 == j5)             e2->at(is, js)->element(j0, j1, j4, j3) +=  2.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j1 == j4)             e2->at(is, js)->element(j5, j2, j0, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j1 == j4 && j2 == j5) e1->at(is, js)->element(j0, j3)         +=  2.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j1 == j5)             e2->at(is, js)->element(j0, j2, j4, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                        if (j1 == j5 && j2 == j4) e1->at(is, js)->element(j0, j3)         += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                      }
                    }
                  }
                }
              }
            }
          }
          for (size_t j7 = 0; j7 != nact; ++j7)
            for (size_t j0 = 0; j0 != nact; ++j0)
              for (size_t j6 = 0; j6 != nact; ++j6)
                for (size_t j5 = 0; j5 != nact; ++j5)
                  for (size_t j2 = 0; j2 != nact; ++j2)
                    for (size_t j1 = 0; j1 != nact; ++j1) {
                      e4->at(is,js)->element(j0, j1, j5, j2, j6, j7) -= 2.0 * Rmat->element(j0, j1, j5, j2, j6, j7);
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        size_t j4i = j4 + nclo;
                        for (size_t j3 = 0; j3 != nact; ++j3) {
                          size_t j3i = j3 + nclo;
                          dshift->element(j3i, j4i) -= Rmat->element(j0, j1, j5, j2, j6, j7) * rdm4->element(j0, j1, j5, j2, j6, j7, j3, j4);
                          if (j4 == j6)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j1, j5, j2, j3, j7);
                          if (j4 == j5)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j1, j3, j2, j6, j7);
                          if (j2 == j6)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j1, j3, j4, j5, j7);
                          if (j2 == j6 && j4 == j5)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j1, j3, j7);
                          if (j2 == j5)                         dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j1, j3, j4, j6, j7);
                          if (j2 == j5 && j4 == j6)             dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j1, j3, j7);
                          if (j2 == j3)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j1, j5, j4, j6, j7);
                          if (j2 == j3 && j4 == j6)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j1, j5, j7);
                          if (j2 == j3 && j4 == j5)             dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j1, j6, j7);
                          if (j1 == j6)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j3, j4, j5, j2, j0, j7);
                          if (j1 == j6 && j4 == j5)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j3, j2, j0, j7);
                          if (j1 == j6 && j2 == j3)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j5, j4, j0, j7);
                          if (j1 == j6 && j2 == j3 && j4 == j5) dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm1->element(j0, j7);
                          if (j1 == j5)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j2, j3, j4, j6, j7);
                          if (j1 == j5 && j4 == j6)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j2, j3, j7);
                          if (j2 == j3 && j1 == j5)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j4, j6, j7);
                          if (j2 == j3 && j1 == j5 && j4 == j6) dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm1->element(j0, j7);
                          if (j1 == j3)                         dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm3->element(j0, j4, j5, j2, j6, j7);
                          if (j1 == j3 && j4 == j6)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j5, j2, j0, j7);
                          if (j1 == j3 && j4 == j5)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j2, j6, j7);
                          if (j2 == j6 && j1 == j3)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j4, j5, j7);
                          if (j2 == j6 && j1 == j3 && j4 == j5) dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm1->element(j0, j7);
                          if (j1 == j3 && j2 == j5)             dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j0, j4, j6, j7);
                          if (j4 == j6 && j2 == j5 && j1 == j3) dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm1->element(j0, j7);
                          if (j1 == j6 && j2 == j5)             dshift->element(j3i, j4i) +=  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j3, j4, j0, j7);
                          if (j1 == j5 && j2 == j6)             dshift->element(j3i, j4i) += -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * rdm2->element(j3, j4, j0, j7);
                          if (j4 == j6)                         e3->at(is,js)->element(j0, j1, j5, j2, j3, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j4 == j5)                         e3->at(is,js)->element(j0, j1, j3, j2, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j6)                         e3->at(is,js)->element(j0, j1, j3, j4, j5, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j6 && j4 == j5)             e2->at(is,js)->element(j0, j1, j3, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j5)                         e3->at(is,js)->element(j0, j1, j3, j4, j6, j7) += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j5 && j4 == j6)             e2->at(is,js)->element(j0, j1, j3, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j3)                         e3->at(is,js)->element(j0, j1, j5, j4, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j3 && j4 == j6)             e2->at(is,js)->element(j0, j1, j5, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j3 && j4 == j5)             e2->at(is,js)->element(j0, j1, j6, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j6)                         e3->at(is,js)->element(j3, j4, j5, j2, j0, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j6 && j4 == j5)             e2->at(is,js)->element(j3, j2, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j6 && j2 == j3)             e2->at(is,js)->element(j5, j4, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j6 && j2 == j3 && j4 == j5) e1->at(is,js)->element(j0, j7)                 += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j5)                         e3->at(is,js)->element(j0, j2, j3, j4, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j5 && j4 == j6)             e2->at(is,js)->element(j0, j2, j3, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j3 && j1 == j5)             e2->at(is,js)->element(j0, j4, j6, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j3 && j1 == j5 && j4 == j6) e1->at(is,js)->element(j0, j7)                 += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j3)                         e3->at(is,js)->element(j0, j4, j5, j2, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j3 && j4 == j6)             e2->at(is,js)->element(j5, j2, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j3 && j4 == j5)             e2->at(is,js)->element(j0, j2, j6, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j6 && j1 == j3)             e2->at(is,js)->element(j0, j4, j5, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j2 == j6 && j1 == j3 && j4 == j5) e1->at(is,js)->element(j0, j7)                 += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j3 && j2 == j5)             e2->at(is,js)->element(j0, j4, j6, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j4 == j6 && j2 == j5 && j1 == j3) e1->at(is,js)->element(j0, j7)                 += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j6 && j2 == j5)             e2->at(is,js)->element(j3, j4, j0, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                          if (j1 == j5 && j2 == j6)             e2->at(is,js)->element(j3, j4, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                        }
                      }
                    }
        }
      }
    }
    ioffset += size_rist;
    timer.tick_print("dshift rist");
  }
    // a i b j
  {
    for (int istate = 0; istate != nstates; ++istate) { // state of T
      shared_ptr<const VectorB> l = lambda[istate];
      shared_ptr<const VectorB> t = amplitude[istate];
      for (int is = 0; is != nstates; ++is) {
        if (info_->sssr() && is != istate) continue;
        double e0loc = e0all_[istate] - e0all_[is];
        for (size_t j3 = 0; j3 != nvirt; ++j3)
          for (size_t j2 = 0; j2 != nclo; ++j2)
            for (size_t j1 = 0; j1 != nvirt; ++j1)
              for (size_t j0 = 0; j0 != nclo; ++j0) {
                const size_t j0i = j0;
                const size_t j1i = j1 + nocc - ncore;
                const size_t j2i = j2;
                const size_t j3i = j3 + nocc - ncore;
                const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3)) + ioffset + (info_->sssr() ? 0 : is * size_aibj);
                const size_t jall2 = j0 + nclo * (j3 + nvirt * (j2 + nclo * j1)) + ioffset + (info_->sssr() ? 0 : is * size_aibj);
                const double lcovar = ((*l)[jall] * 8.0 - (*l)[jall2] * 4.0);
                const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j1+nocc] + eig_[j3+nocc] + e0loc;
                const double Lambda = shift2 * lcovar * (*t)[jall];
                dshift->element(j0i, j0i) += Lambda;
                dshift->element(j1i, j1i) -= Lambda;
                dshift->element(j2i, j2i) += Lambda;
                dshift->element(j3i, j3i) -= Lambda;
                nimag[istate] += Lambda;
                nimag[is] -= Lambda;
                ((*e0->at(istate,istate))) += 2.0 * lcovar * (*t)[jall] * denom;
              }
      }
    }
    ioffset += size_aibj * (info_->sssr() ? 1 : nstates);
    timer.tick_print("dshift aibj");
  }

  return tie(dshift, e0, e1, e2, e3, e4, nimag);
}

#endif
