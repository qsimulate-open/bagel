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


tuple<shared_ptr<Matrix>,shared_ptr<VecRDM<1>>,shared_ptr<VecRDM<2>>,shared_ptr<VecRDM<3>>,shared_ptr<VecRDM<3>>,vector<double>>
  MSCASPT2::MSCASPT2::make_d2_imag() const {
  const size_t nact = info_->nact();
  const int nstates = nact ? info_->ciwfn()->nstates() : 1;
  const size_t nclosed = info_->nclosed();
  const size_t ncore = info_->ncore();
  const size_t nclo = nclosed - ncore;

  Timer timer(1);

  shared_ptr<Matrix> dshift = den2_->clone();
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
          auto e1temp = make_shared<RDM<1>>(nact);
          auto e2temp = make_shared<RDM<2>>(nact);
          auto e3temp = make_shared<RDM<3>>(nact);
          auto e4temp = make_shared<RDM<3>>(nact);

          e1->emplace(is, js, e1temp);
          e2->emplace(is, js, e2temp);
          e3->emplace(is, js, e3temp);
          e4->emplace(is, js, e4temp);
        }
  }

  const double shift2 = info_->shift() * info_->shift();

  for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
    auto smallzall = make_shared<Vec<Matrix>>();
    auto largexall = make_shared<Vec<Matrix>>();
    auto largeyall = make_shared<Vec<Matrix>>();
    auto largeqall = make_shared<Vec<Matrix>>();

    for (int istate = 0; istate != nstates; ++istate) {
      if (!info_->sssr() && istate != 0) continue;        // if ms-mr, only one needed
      const size_t interm_size = (iext != Excitations::aibj) ? t_orthogonal_->shalf(iext, istate).ndim() : 1;
      auto temp = make_shared<Matrix>(interm_size, interm_size);
      smallzall->emplace(istate, temp->copy());
      largexall->emplace(istate, temp->copy());
      largeyall->emplace(istate, temp->copy());
      largeqall->emplace(istate, temp->copy());
    }

    switch (iext) {
      case Excitations::arbs:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i3 : t_orthogonal_->virt())
            for (auto& i1 : t_orthogonal_->virt())
              for (auto& i0o : t_orthogonal_->interm(dataindex)) {
                if (!t->is_local(i0o, i1, i3)) continue;
                const unique_ptr<double[]> amplitude = t->get_block(i0o, i1, i3);
                const unique_ptr<double[]> lambda    = l->get_block(i0o, i1, i3);
                const unique_ptr<double[]> denom     = d->get_block(i0o, i1, i3);
                for (size_t j3 = i3.offset(), jall = 0; j3 != i3.offset()+i3.size(); ++j3) {
                  const size_t j3i = j3 - ncore;
                  for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                    const size_t j1i = j1 - ncore;
                    for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                      const double Lambda = shift2 * lambda[jall] * amplitude[jall];
                      dshift->element(j1i, j1i) -= Lambda;
                      dshift->element(j3i, j3i) -= Lambda;
                      smallz->element(j0o, j0o) -= Lambda;
                      nimag[istate] -= Lambda;
                    }
                  }
                }
                for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                  const unique_ptr<double[]> amplitudek = t->get_block(i1o, i1, i3);
                  const unique_ptr<double[]> lambdak    = l->get_block(i1o, i1, i3);
                  const unique_ptr<double[]> denomk     = d->get_block(i1o, i1, i3);
                  for (size_t j3 = i3.offset(), jall = 0; j3 != i3.offset()+i3.size(); ++j3) {
                    for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                      for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                        for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                          const size_t kall = j1o - i1o.offset() + i1o.size() * (j1 - i1.offset() + i1.size() * (j3 - i3.offset()));
                          const double lt = lambda[jall] * amplitudek[kall] * shift2 * denom[jall];
                          const double tl = amplitude[jall] * lambdak[kall] * shift2 * denomk[kall];
                          largey->element(j0o, j1o) += lambda[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                          largeq->element(j0o, j1o) += (lt + tl);
                        }
                      }
                    }
                  }
                }
              }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }

        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<2>>(nact);
              // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                        const double VrsO = VmatJ.element(j0o, j2 + j3 * nact);
                        Rmat->element(j0, j1, j2, j3) += Qmat->element(j0o, j0 + j1 * nact) * VrsO;
                        e2->at(js, is)->element(j0, j2, j1, j3) += Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                      }
                    }
                  }
                }
              }
              Rmat->allreduce();
              // (3) form d^{(2)}_{rs} and e3
              for (size_t j0 = 0; j0 != nact; ++j0)
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2)
                    for (size_t j3 = 0; j3 != nact; ++j3)
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        const size_t j4i = j4 + nclo;
                        for (size_t j5 = 0; j5 != nact; ++j5) {
                          const size_t j5i = j5 + nclo;
                          dshift->element(j4i, j5i) += Rmat->element(j0, j1, j2, j3) * rdm3->element(j0, j2, j1, j3, j4, j5);
                          e3->at(js, is)->element(j0, j2, j1, j3, j4, j5) += 2.0 * Rmat->element(j0, j1, j2, j3) * fockact_->element(j4, j5);
                        }
                      }
                }
            }
          }
        }
        timer.tick_print("dshift arbs");
        break;
      case Excitations::arbi:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i3 : t_orthogonal_->virt())
            for (auto& i2 : t_orthogonal_->closed())
              for (auto& i1 : t_orthogonal_->virt())
                for (auto& i0o : t_orthogonal_->interm(dataindex)) {
                  if (!t->is_local(i0o, i1, i2, i3)) continue;
                  const unique_ptr<double[]> amplitude = t->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> lambda    = l->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> denom     = d->get_block(i0o, i1, i2, i3);
                  unique_ptr<double[]> tcovar          = t->get_block(i0o, i1, i2, i3);
                  unique_ptr<double[]> lcovar          = l->get_block(i0o, i1, i2, i3);
                  {
                    const unique_ptr<double[]> amplitude2 = t->get_block(i0o, i3, i2, i1);
                    const unique_ptr<double[]> lambda2    = l->get_block(i0o, i3, i2, i1);
                    sort_indices<0,3,2,1,2,1,-1,1>(amplitude2, tcovar, i0o.size(), i3.size(), i2.size(), i1.size());
                    sort_indices<0,3,2,1,2,1,-1,1>(lambda2, lcovar, i0o.size(), i3.size(), i2.size(), i1.size());
                  }
                  for (size_t j3 = i3.offset(), jall = 0; j3 != i3.offset()+i3.size(); ++j3) {
                    const size_t j3i = j3 - ncore;
                    for (size_t j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) {
                      const size_t j2i = j2 - ncore;
                      for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                        const size_t j1i = j1 - ncore;
                        for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                          const double Lambda = shift2 * lcovar[jall] * amplitude[jall];
                          dshift->element(j1i, j1i) -= Lambda;
                          dshift->element(j2i, j2i) += Lambda;
                          dshift->element(j3i, j3i) -= Lambda;
                          smallz->element(j0o, j0o) -= Lambda;
                          nimag[istate] -= Lambda;
                        }
                      }
                    }
                  }
                  for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                    const unique_ptr<double[]> amplitudek = t->get_block(i1o, i1, i2, i3);
                    const unique_ptr<double[]> lambdak    = l->get_block(i1o, i1, i2, i3);
                    const unique_ptr<double[]> denomk     = d->get_block(i1o, i1, i2, i3);
                    for (size_t j3 = i3.offset(), jall = 0; j3 != i3.offset()+i3.size(); ++j3) {
                      for (size_t j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) {
                        for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                          for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                            for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                              const size_t kall = j1o - i1o.offset() + i1o.size() * (j1 - i1.offset() + i1.size() * (j2 - i2.offset() + i2.size() * (j3 - i3.offset())));
                              const double lt = lcovar[jall] * amplitudek[kall] * shift2 * denom[jall];
                              const double tl = tcovar[jall] * lambdak[kall] * shift2 * denomk[kall];
                              largey->element(j0o, j1o) += lcovar[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                              largeq->element(j0o, j1o) += (lt + tl);
                            }
                          }
                        }
                      }
                    }
                  }
                }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<1>>(nact);
              // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                    const double VrO = VmatJ.element(j0o, j1);
                    Rmat->element(j0, j1) += Qmat->element(j0o, j0) * VrO;
                    e1->at(js, is)->element(j0, j1) += Pmat->element(j0o, j0) * VrO;
                  }
                }
              }
              Rmat->allreduce();
              // (3) form d^{(2)}_{rs} and e3
              for (size_t j0 = 0; j0 != nact; ++j0)
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    const size_t j2i = j2 + nclo;
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      const size_t j3i = j3 + nclo;
                      dshift->element(j2i, j3i) += Rmat->element(j0, j1) * rdm2->element(j0, j1, j2, j3);
                      e2->at(js, is)->element(j0, j1, j2, j3) += 2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3);
                    }
                  }
                }
            }
          }
        }
        timer.tick_print("dshift arbi");
        break;
      case Excitations::airj:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i2 : t_orthogonal_->closed())
            for (auto& i1 : t_orthogonal_->virt())
              for (auto& i0 : t_orthogonal_->closed())
                for (auto& i0o : t_orthogonal_->interm(dataindex)) {
                  if (!t->is_local(i0o, i0, i1, i2)) continue;
                  const unique_ptr<double[]> amplitude = t->get_block(i0o, i0, i1, i2);
                  const unique_ptr<double[]> lambda    = l->get_block(i0o, i0, i1, i2);
                  const unique_ptr<double[]> denom     = d->get_block(i0o, i0, i1, i2);
                  unique_ptr<double[]> tcovar          = t->get_block(i0o, i0, i1, i2);
                  unique_ptr<double[]> lcovar          = l->get_block(i0o, i0, i1, i2);
                  {
                    const unique_ptr<double[]> amplitude2 = t->get_block(i0o, i2, i1, i0);
                    const unique_ptr<double[]> lambda2    = l->get_block(i0o, i2, i1, i0);
                    sort_indices<0,3,2,1,2,1,-1,1>(amplitude2, tcovar, i0o.size(), i2.size(), i1.size(), i0.size());
                    sort_indices<0,3,2,1,2,1,-1,1>(lambda2, lcovar, i0o.size(), i2.size(), i1.size(), i0.size());
                  }
                  for (size_t j2 = i2.offset(), jall = 0; j2 != i2.offset()+i2.size(); ++j2) {
                    const size_t j2i = j2 - ncore;
                    for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                      const size_t j1i = j1 - ncore;
                      for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                        const size_t j0i = j0 - ncore;
                        for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                          const double Lambda = shift2 * lcovar[jall] * amplitude[jall];
                          dshift->element(j0i, j0i) += Lambda;
                          dshift->element(j1i, j1i) -= Lambda;
                          dshift->element(j2i, j2i) += Lambda;
                          smallz->element(j0o, j0o) -= Lambda;
                          nimag[istate] -= Lambda;
                        }
                      }
                    }
                  }
                  for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                    const unique_ptr<double[]> amplitudek = t->get_block(i1o, i0, i1, i2);
                    const unique_ptr<double[]> lambdak    = l->get_block(i1o, i0, i1, i2);
                    const unique_ptr<double[]> denomk     = d->get_block(i1o, i0, i1, i2);
                    for (size_t j2 = i2.offset(), jall = 0; j2 != i2.offset()+i2.size(); ++j2) {
                      for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                        for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                          for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                            for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                              const size_t kall = j1o - i1o.offset() + i1o.size() * (j0 - i0.offset() + i0.size() * (j1 - i1.offset() + i1.size() * (j2 - i2.offset())));
                              const double lt = lcovar[jall] * amplitudek[kall] * shift2 * denom[jall];
                              const double tl = tcovar[jall] * lambdak[kall] * shift2 * denomk[kall];
                              largey->element(j0o, j1o) += lcovar[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                              largeq->element(j0o, j1o) += (lt + tl);
                            }
                          }
                        }
                      }
                    }
                  }
                }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<1>>(nact);
              // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                    const double VrO = VmatJ.element(j0o, j1);
                    Rmat->element(j0, j1) += Qmat->element(j0o, j0) * VrO;
                    e1->at(js, is)->element(j0, j1) += -Pmat->element(j0o, j0) * VrO;
                  }
                }
              }
              Rmat->allreduce();
              // (3) form d^{(2)}_{rs} and e3
              for (size_t j0 = 0; j0 != nact; ++j0)
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    const size_t j2i = j2 + nclo;
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      const size_t j3i = j3 + nclo;
                      dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm2->element(j0, j1, j2, j3);
                      if (j1 == j2 && j0 == j3 && is == js) dshift->element(j2i, j3i) +=  2.0 * Rmat->element(j0, j1);
                      if (j1 == j2)                         dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm1->element(j0, j3);
                      if (j0 == j3)                         dshift->element(j2i, j3i) += -1.0 * Rmat->element(j0, j1) * rdm1->element(j1, j2);
                      if (j0 == j1 && is == js)             dshift->element(j2i, j3i) +=  2.0 * Rmat->element(j0, j1) * rdm1->element(j2, j3);
                      e2->at(js, is)->element(j0, j1, j2, j3) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                      if (j1 == j2)                         e1->at(js, is)->element(j0, j3) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                      if (j0 == j3)                         e1->at(js, is)->element(j1, j2) += -1.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                      if (j0 == j1 && is == js)             e1->at(js, is)->element(j2, j3) +=  2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3) * 2.0;
                    }
                  }
                }
            }
          }
        }
        timer.tick_print("dshift airj");
        break;
      case Excitations::risj:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i2 : t_orthogonal_->closed())
            for (auto& i0 : t_orthogonal_->closed())
              for (auto& i0o : t_orthogonal_->interm(dataindex)) {
                if (!t->is_local(i0o, i0, i2)) continue;
                const unique_ptr<double[]> amplitude = t->get_block(i0o, i0, i2);
                const unique_ptr<double[]> lambda    = l->get_block(i0o, i0, i2);
                const unique_ptr<double[]> denom     = d->get_block(i0o, i0, i2);
                for (size_t j2 = i2.offset(), jall = 0; j2 != i2.offset()+i2.size(); ++j2) {
                  const size_t j2i = j2 - ncore;
                  for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                    const size_t j0i = j0 - ncore;
                    for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                      const double Lambda = shift2 * lambda[jall] * amplitude[jall];
                      dshift->element(j0i, j0i) += Lambda;
                      dshift->element(j2i, j2i) += Lambda;
                      smallz->element(j0o, j0o) -= Lambda;
                      nimag[istate] -= Lambda;
                    }
                  }
                }
                for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                  const unique_ptr<double[]> amplitudek = t->get_block(i1o, i0, i2);
                  const unique_ptr<double[]> lambdak    = l->get_block(i1o, i0, i2);
                  const unique_ptr<double[]> denomk     = d->get_block(i1o, i0, i2);
                  for (size_t j2 = i2.offset(), jall = 0; j2 != i2.offset()+i2.size(); ++j2) {
                    for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                      for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                        for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                          const size_t kall = j1o - i1o.offset() + i1o.size() * (j0 - i0.offset() + i0.size() * (j2 - i2.offset()));
                          const double lt = lambda[jall] * amplitudek[kall] * shift2 * denom[jall];
                          const double tl = amplitude[jall] * lambdak[kall] * shift2 * denomk[kall];
                          largey->element(j0o, j1o) += lambda[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                          largeq->element(j0o, j1o) += (lt + tl);
                        }
                      }
                    }
                  }
                }
              }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<2>>(nact);
              // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                        const double VrsO = VmatJ.element(j0o, j2 + j3 * nact);
                        Rmat->element(j0, j2, j1, j3) += Qmat->element(j0o, j0 + j1 * nact) * VrsO;
                        e2->at(js, is)->element(j0, j2, j1, j3) += Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                        if (j2 == j1) e1->at(js, is)->element(j0, j3) +=  1.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                        if (j1 == j3) e1->at(js, is)->element(j0, j2) += -2.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                        if (j0 == j2) e1->at(js, is)->element(j1, j3) += -2.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                        if (j0 == j3) e1->at(js, is)->element(j1, j2) +=  1.0 * Pmat->element(j0o, j0 + j1 * nact) * VrsO;
                      }
                    }
                  }
                }
              }
              Rmat->allreduce();
              // (3) form d^{(2)}_{rs} and e3
              for (size_t j4 = 0; j4 != nact; ++j4)
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j4;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
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
                          e3->at(js,is)->element(j0, j5, j1, j4, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j2 == j5)                                     e2->at(js,is)->element(j1, j4, j0, j3) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j2 == j4)                                     e2->at(js,is)->element(j0, j5, j1, j3) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j5)                                     e2->at(js,is)->element(j0, j4, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j5 && j2 == j4)                         e1->at(js,is)->element(j0, j3)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j4)                                     e2->at(js,is)->element(j0, j5, j2, j3) += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j4 && j2 == j5)                         e1->at(js,is)->element(j0, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j3)                                     e2->at(js,is)->element(j0, j5, j2, j4) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j3 && j2 == j5)                         e1->at(js,is)->element(j0, j4)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j3 && j2 == j4)                         e1->at(js,is)->element(j0, j5)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j5)                                     e2->at(js,is)->element(j1, j4, j2, j3) += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j5 && j2 == j4)                         e1->at(js,is)->element(j1, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j3 && j0 == j5)                         e1->at(js,is)->element(j2, j4)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j4)                                     e2->at(js,is)->element(j1, j5, j2, j3) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j4 && j2 == j5)                         e1->at(js,is)->element(j1, j3)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j3 && j0 == j4)                         e1->at(js,is)->element(j2, j5)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j3)                                     e2->at(js,is)->element(j2, j5, j1, j4) += 2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j3 && j2 == j5)                         e1->at(js,is)->element(j1, j4)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j3 && j2 == j4)                         e1->at(js,is)->element(j1, j5)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j1 == j5 && j0 == j3)                         e1->at(js,is)->element(j2, j4)         += 2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j3 && j1 == j4)                         e1->at(js,is)->element(j2, j5)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j5 && j1 == j4)                         e1->at(js,is)->element(j2, j3)         += 2.0 *  4.0 * factor * fockact_->element(j2, j3);
                          if (j0 == j4 && j1 == j5)                         e1->at(js,is)->element(j2, j3)         += 2.0 * -2.0 * factor * fockact_->element(j2, j3);
                        }
                      }
                }
            }
          }
        }
        timer.tick_print("dshift risj");
        break;
      case Excitations::airs:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i1 : t_orthogonal_->virt())
            for (auto& i0 : t_orthogonal_->closed())
              for (auto& i0o : t_orthogonal_->interm(dataindex)) {
                if (!t->is_local(i0o, i0, i1)) continue;
                const unique_ptr<double[]> amplitude = t->get_block(i0o, i0, i1);
                const unique_ptr<double[]> lambda    = l->get_block(i0o, i0, i1);
                const unique_ptr<double[]> denom     = d->get_block(i0o, i0, i1);
                for (size_t j1 = i1.offset(), jall = 0; j1 != i1.offset()+i1.size(); ++j1) {
                  const size_t j1i = j1 - ncore;
                  for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                    const size_t j0i = j0 - ncore;
                    for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                      const double Lambda = shift2 * lambda[jall] * amplitude[jall];
                      dshift->element(j0i, j0i) += Lambda;
                      dshift->element(j1i, j1i) -= Lambda;
                      smallz->element(j0o, j0o) -= Lambda;
                      nimag[istate] -= Lambda;
                    }
                  }
                }
                for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                  const unique_ptr<double[]> amplitudek = t->get_block(i1o, i0, i1);
                  const unique_ptr<double[]> lambdak    = l->get_block(i1o, i0, i1);
                  const unique_ptr<double[]> denomk     = d->get_block(i1o, i0, i1);
                  for (size_t j1 = i1.offset(), jall = 0; j1 != i1.offset()+i1.size(); ++j1) {
                    for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) {
                      for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                        for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                          const size_t kall = j1o - i1o.offset() + i1o.size() * (j0 - i0.offset() + i0.size() * (j1 - i1.offset()));
                          const double lt = lambda[jall] * amplitudek[kall] * shift2 * denom[jall];
                          const double tl = amplitude[jall] * lambdak[kall] * shift2 * denomk[kall];
                          largey->element(j0o, j1o) += lambda[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                          largeq->element(j0o, j1o) += (lt + tl);
                        }
                      }
                    }
                  }
                }
              }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            shared_ptr<RDM<1>> rdm1c;
            shared_ptr<RDM<2>> rdm2c;
            shared_ptr<RDM<3>> rdm3c;
            shared_ptr<RDM<4>> rdm4c;
            tie(rdm1c, rdm2c, rdm3c, rdm4c) = feed_rdm(is, js);
            {
              shared_ptr<Matrix> VmatAI, VmatBI;
              shared_ptr<Matrix> VmatAJ, VmatBJ;
              {
                auto VmatI = make_shared<Matrix>(t_orthogonal_->shalf(iext, is));
                VmatAI = VmatI->get_submatrix(0,           0, interm_size, nact * nact);
                VmatBI = VmatI->get_submatrix(0, nact * nact, interm_size, nact * nact);
                auto VmatJ = make_shared<Matrix>(t_orthogonal_->shalf(iext, js));
                VmatAJ = VmatJ->get_submatrix(0,           0, interm_size, nact * nact);
                VmatBJ = VmatJ->get_submatrix(0, nact * nact, interm_size, nact * nact);
              }
              auto QmatA = make_shared<Matrix>(*smallz % (2.0 * (*VmatAI) - (*VmatBI)));
              auto QmatB = make_shared<Matrix>(*smallz % (-1.0 * (*VmatAI)));
              auto QmatC = make_shared<Matrix>(*smallz % (*VmatBI));
              auto PmatA = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (2.0 * (*VmatAI) - (*VmatBI)));
              auto PmatB = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (-1.0 * (*VmatAI)));
              auto PmatC = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*VmatBI));
              auto RmatA = make_shared<RDM<2>>(nact);
              auto RmatC = make_shared<RDM<2>>(nact);
              // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                        const double VrsO = VmatAJ->element(j0o, j2 + j3 * nact);
                        const double VrsS = VmatBJ->element(j0o, j2 + j3 * nact);
                        RmatA->element(j0, j1, j3, j2) += QmatA->element(j0o, j0 + j1 * nact) * VrsO + QmatB->element(j0o, j0 + j1 * nact) * VrsS;
                        RmatC->element(j0, j1, j3, j2) += QmatC->element(j0o, j0 + j1 * nact) * VrsS;
                        e2->at(js, is)->element(j0, j1, j3, j2)       += (PmatA->element(j0o, j0 + j1 * nact) * VrsO + PmatB->element(j0o, j0 + j1 * nact) * VrsS);
                        if (j1 == j3) e1->at(js, is)->element(j0, j2) += (PmatA->element(j0o, j0 + j1 * nact) * VrsO + PmatB->element(j0o, j0 + j1 * nact) * VrsS);
                        e2->at(is, js)->element(j1, j3, j2, j0)       += -1.0 * PmatC->element(j0o, j0 + j1 * nact) * VrsS;
                        if (j1 == j3) e1->at(is, js)->element(j2, j0) +=  2.0 * PmatC->element(j0o, j0 + j1 * nact) * VrsS;
                      }
                    }
                  }
                }
              }
              RmatA->allreduce();
              RmatC->allreduce();
              // (3) form d^{(2)}_{rs} and e3
              for (size_t j4 = 0; j4 != nact; ++j4)
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j4;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
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
                          dshift->element(j2i, j3i) -= RmatC->element(j0, j1, j4, j5) * rdm3c->element(j1, j4, j5, j0, j2, j3);
                          if (j3 == j1)             dshift->element(j2i, j3i) += -1.0 * RmatC->element(j0, j1, j4, j5) * rdm2c->element(j2, j4, j5, j0);
                          if (j4 == j2)             dshift->element(j2i, j3i) += -1.0 * RmatC->element(j0, j1, j4, j5) * rdm2c->element(j1, j3, j5, j0);
                          if (j3 == j1 && j4 == j2) dshift->element(j2i, j3i) +=  2.0 * RmatC->element(j0, j1, j4, j5) * rdm1c->element(j5, j0);
                          if (j4 == j1)             dshift->element(j2i, j3i) +=  2.0 * RmatC->element(j0, j1, j4, j5) * rdm2c->element(j2, j3, j5, j0);
                          e3->at(js,is)->element(j0, j1, j4, j5, j2, j3) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                          if (j3 == j4)             e2->at(js,is)->element(j0, j1, j2, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                          if (j1 == j2)             e2->at(js,is)->element(j0, j3, j4, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                          if (j1 == j2 && j3 == j4) e1->at(js,is)->element(j0, j5)         += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
                          if (j1 == j4)             e2->at(js,is)->element(j2, j3, j0, j5) += 2.0 * RmatA->element(j0, j1, j4, j5) * fockact_->element(j2, j3);
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
        }
        timer.tick_print("dshift airs");
        break;
      case Excitations::arst:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i0 : t_orthogonal_->virt())
            for (auto& i0o : t_orthogonal_->interm(dataindex)) {
              if (!t->is_local(i0o, i0)) continue;
              const unique_ptr<double[]> amplitude = t->get_block(i0o, i0);
              const unique_ptr<double[]> lambda    = l->get_block(i0o, i0);
              const unique_ptr<double[]> denom     = d->get_block(i0o, i0);
              for (size_t j0 = i0.offset(), jall = 0; j0 != i0.offset()+i0.size(); ++j0) {
                const size_t j0i = j0 - ncore;
                for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                  const double Lambda = shift2 * lambda[jall] * amplitude[jall];
                  dshift->element(j0i, j0i) -= Lambda;
                  smallz->element(j0o, j0o) -= Lambda;
                  nimag[istate] -= Lambda;
                }
              }
              for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                const unique_ptr<double[]> amplitudek = t->get_block(i1o, i0);
                const unique_ptr<double[]> lambdak    = l->get_block(i1o, i0);
                const unique_ptr<double[]> denomk     = d->get_block(i1o, i0);
                for (size_t j0 = i0.offset(), jall = 0; j0 != i0.offset()+i0.size(); ++j0) {
                  for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                    for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                      const size_t kall = j1o - i1o.offset() + i1o.size() * (j0 - i0.offset());
                      const double lt = lambda[jall] * amplitudek[kall] * shift2 * denom[jall];
                      const double tl = amplitude[jall] * lambdak[kall] * shift2 * denomk[kall];
                      largey->element(j0o, j1o) += lambda[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                      largeq->element(j0o, j1o) += (lt + tl);
                    }
                  }
                }
              }
            }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<3>>(nact);
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j1 + nact * j0;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        for (size_t j5 = 0; j5 != nact; ++j5) {
                          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                            const double VrstO = VmatJ.element(j0o, j3 + nact * (j4 + nact * j5));
                            Rmat->element(j1, j2, j5, j4, j0, j3) += Qmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            e3->at(js, is)->element(j0, j1, j5, j4, j2, j3) += Pmat->element(j0o, j2 + nact * (j0 + nact * j1)) * VrstO;
                            if (j1 == j5) e2->at(js, is)->element(j0, j4, j2, j3) += Pmat->element(j0o, j2 + nact * (j0 + nact * j1)) * VrstO;
                          }
                        }
                      }
                    }
                  }
                }
              }
              Rmat->allreduce();
              for (size_t j7 = 0; j7 != nact; ++j7)
                for (size_t j0 = 0; j0 != nact; ++j0) {
                  const int jall = j0 + nact * j7;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j6 = 0; j6 != nact; ++j6)
                    for (size_t j5 = 0; j5 != nact; ++j5)
                      for (size_t j2 = 0; j2 != nact; ++j2)
                        for (size_t j1 = 0; j1 != nact; ++j1) {
                          e4->at(js,is)->element(j1, j2, j5, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7);
                          for (size_t j4 = 0; j4 != nact; ++j4) {
                            size_t j4i = j4 + nclo;
                            for (size_t j3 = 0; j3 != nact; ++j3) {
                              size_t j3i = j3 + nclo;
                              dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm4->element(j1, j2, j5, j6, j0, j7, j3, j4);
                              if (j4 == j5)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j1, j2, j3, j6, j0, j7);
                              if (j2 == j3)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j1, j4, j5, j6, j0, j7);
                              if (j2 == j3 && j4 == j5) dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm2->element(j1, j6, j0, j7);
                              if (j2 == j5)             dshift->element(j3i, j4i) += Rmat->element(j1, j2, j5, j6, j0, j7) * rdm3->element(j3, j4, j1, j6, j0, j7);
                              if (j4 == j5)             e3->at(js,is)->element(j1, j2, j3, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                              if (j2 == j3)             e3->at(js,is)->element(j1, j4, j5, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                              if (j2 == j3 && j4 == j5) e2->at(js,is)->element(j1, j6, j0, j7)         += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                              if (j2 == j5)             e3->at(js,is)->element(j3, j4, j1, j6, j0, j7) += 2.0 * Rmat->element(j1, j2, j5, j6, j0, j7) * fockact_->element(j3, j4);
                            }
                          }
                        }
                }
            }
          }
        }
        timer.tick_print("dshift arst");
        break;
      case Excitations::rist:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(iext);
          const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(iext);

          const size_t dataindex = info_->sssr() ? iext + istate * Excitations::total : iext;
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(istate) : largeqall->at(0);

          for (auto& i0 : t_orthogonal_->closed())
            for (auto& i0o : t_orthogonal_->interm(dataindex)) {
              if (!t->is_local(i0o, i0)) continue;
              const unique_ptr<double[]> amplitude = t->get_block(i0o, i0);
              const unique_ptr<double[]> lambda    = l->get_block(i0o, i0);
              const unique_ptr<double[]> denom     = d->get_block(i0o, i0);
              for (size_t j0 = i0.offset(), jall = 0; j0 != i0.offset()+i0.size(); ++j0) {
                const size_t j0i = j0 - ncore;
                for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                  const double Lambda = shift2 * lambda[jall] * amplitude[jall];
                  dshift->element(j0i, j0i) += Lambda;
                  smallz->element(j0o, j0o) -= Lambda;
                  nimag[istate] -= Lambda;
                }
              }
              for (auto& i1o : t_orthogonal_->interm(dataindex)) {
                const unique_ptr<double[]> amplitudek = t->get_block(i1o, i0);
                const unique_ptr<double[]> lambdak    = l->get_block(i1o, i0);
                const unique_ptr<double[]> denomk     = d->get_block(i1o, i0);
                for (size_t j0 = i0.offset(), jall = 0; j0 != i0.offset()+i0.size(); ++j0) {
                  for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++jall) {
                    for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                      const size_t kall = j1o - i1o.offset() + i1o.size() * (j0 - i0.offset());
                      const double lt = lambda[jall] * amplitudek[kall] * shift2 * denom[jall];
                      const double tl = amplitude[jall] * lambdak[kall] * shift2 * denomk[kall];
                      largey->element(j0o, j1o) += lambda[jall] * amplitudek[kall] * shift2 * (denom[jall] - denomk[kall]);
                      largeq->element(j0o, j1o) += (lt + tl);
                    }
                  }
                }
              }
            }
        }
        for (int istate = 0; istate != nstates; ++istate) {
          if (!info_->sssr() && istate != 0) continue;
          smallzall->at(istate)->allreduce();
          largeyall->at(istate)->allreduce();
          largeqall->at(istate)->allreduce();

          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(istate) : smallzall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(istate) : largexall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(istate) : largeyall->at(0);

          const size_t interm_size = t_orthogonal_->shalf(iext, istate).ndim();

          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              if (j1o == j0o) continue;
              const double fdiff = t_orthogonal_->phi(iext, istate, j1o) - t_orthogonal_->phi(iext, istate, j0o);
              smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
            }
          }
          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
            for (size_t j1o = 0; j1o != interm_size; ++j1o) {
              largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * t_orthogonal_->phi(iext, istate, j1o))
                                        + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * t_orthogonal_->phi(iext, istate, j0o));
            }
          }
        }
        for (size_t is = 0; is != nstates; ++is) {
          shared_ptr<Matrix> smallz = info_->sssr() ? smallzall->at(is) : smallzall->at(0);
          shared_ptr<Matrix> largey = info_->sssr() ? largeyall->at(is) : largeyall->at(0);
          shared_ptr<Matrix> largeq = info_->sssr() ? largeqall->at(is) : largeqall->at(0);
          shared_ptr<Matrix> largex = info_->sssr() ? largexall->at(is) : largexall->at(0);
          const size_t interm_size = t_orthogonal_->shalf(iext, is).ndim();
          for (size_t js = 0; js != nstates; ++js) {
            if (info_->sssr() && is != js) continue;
            shared_ptr<RDM<1>> rdm1;
            shared_ptr<RDM<2>> rdm2;
            shared_ptr<RDM<3>> rdm3;
            shared_ptr<RDM<4>> rdm4;
            tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
            {
              MatView VmatI = t_orthogonal_->shalf(iext, is);
              MatView VmatJ = t_orthogonal_->shalf(iext, js);
              auto Qmat = make_shared<Matrix>((*smallz) % VmatI);
              auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % VmatI);
              auto Rmat = make_shared<RDM<3>>(nact);
              for (size_t j0 = 0; j0 != nact; ++j0) {
                for (size_t j1 = 0; j1 != nact; ++j1) {
                  const int jall = j0 + nact * j1;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j2 = 0; j2 != nact; ++j2) {
                    for (size_t j3 = 0; j3 != nact; ++j3) {
                      for (size_t j4 = 0; j4 != nact; ++j4) {
                        for (size_t j5 = 0; j5 != nact; ++j5) {
                          for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                            const double VrstO = VmatJ.element(j0o, j3 + nact * (j4 + nact * j5));
                            Rmat->element(j0, j1, j5, j2, j4, j3) += Qmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            e3->at(js, is)->element(j0, j1, j5, j2, j4, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j2 == j4)             e2->at(js, is)->element(j0, j1, j5, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j2 == j5)             e2->at(js, is)->element(j0, j1, j4, j3) +=  2.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j1 == j4)             e2->at(js, is)->element(j5, j2, j0, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j1 == j4 && j2 == j5) e1->at(js, is)->element(j0, j3)         +=  2.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j1 == j5)             e2->at(js, is)->element(j0, j2, j4, j3) += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                            if (j1 == j5 && j2 == j4) e1->at(js, is)->element(j0, j3)         += -1.0 * Pmat->element(j0o, j0 + nact * (j1 + nact * j2)) * VrstO;
                          }
                        }
                      }
                    }
                  }
                }
              }
              Rmat->allreduce();
              for (size_t j7 = 0; j7 != nact; ++j7)
                for (size_t j0 = 0; j0 != nact; ++j0) {
                  const int jall = j0 + nact * j7;
                  if (jall % mpi__->size() != mpi__->rank()) continue;
                  for (size_t j6 = 0; j6 != nact; ++j6)
                    for (size_t j5 = 0; j5 != nact; ++j5)
                      for (size_t j2 = 0; j2 != nact; ++j2)
                        for (size_t j1 = 0; j1 != nact; ++j1) {
                          e4->at(js,is)->element(j0, j1, j5, j2, j6, j7) -= 2.0 * Rmat->element(j0, j1, j5, j2, j6, j7);
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
                              if (j4 == j6)                         e3->at(js,is)->element(j0, j1, j5, j2, j3, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j4 == j5)                         e3->at(js,is)->element(j0, j1, j3, j2, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j6)                         e3->at(js,is)->element(j0, j1, j3, j4, j5, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j6 && j4 == j5)             e2->at(js,is)->element(j0, j1, j3, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j5)                         e3->at(js,is)->element(j0, j1, j3, j4, j6, j7) += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j5 && j4 == j6)             e2->at(js,is)->element(j0, j1, j3, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j3)                         e3->at(js,is)->element(j0, j1, j5, j4, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j3 && j4 == j6)             e2->at(js,is)->element(j0, j1, j5, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j3 && j4 == j5)             e2->at(js,is)->element(j0, j1, j6, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j6)                         e3->at(js,is)->element(j3, j4, j5, j2, j0, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j6 && j4 == j5)             e2->at(js,is)->element(j3, j2, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j6 && j2 == j3)             e2->at(js,is)->element(j5, j4, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j6 && j2 == j3 && j4 == j5) e1->at(js,is)->element(j0, j7)                 += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j5)                         e3->at(js,is)->element(j0, j2, j3, j4, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j5 && j4 == j6)             e2->at(js,is)->element(j0, j2, j3, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j3 && j1 == j5)             e2->at(js,is)->element(j0, j4, j6, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j3 && j1 == j5 && j4 == j6) e1->at(js,is)->element(j0, j7)                 += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j3)                         e3->at(js,is)->element(j0, j4, j5, j2, j6, j7) += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j3 && j4 == j6)             e2->at(js,is)->element(j5, j2, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j3 && j4 == j5)             e2->at(js,is)->element(j0, j2, j6, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j6 && j1 == j3)             e2->at(js,is)->element(j0, j4, j5, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j2 == j6 && j1 == j3 && j4 == j5) e1->at(js,is)->element(j0, j7)                 += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j3 && j2 == j5)             e2->at(js,is)->element(j0, j4, j6, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j4 == j6 && j2 == j5 && j1 == j3) e1->at(js,is)->element(j0, j7)                 += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j6 && j2 == j5)             e2->at(js,is)->element(j3, j4, j0, j7)         += 2.0 *  2.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                              if (j1 == j5 && j2 == j6)             e2->at(js,is)->element(j3, j4, j0, j7)         += 2.0 * -1.0 * Rmat->element(j0, j1, j5, j2, j6, j7) * fockact_->element(j3, j4);
                            }
                          }
                        }
                }
            }
          }
        }
        timer.tick_print("dshift rist");
        break;
      case Excitations::aibj:
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          for (int is = 0; is != nstates; ++is) {
            if (info_->sssr() && is != istate) continue;
            const int pos = iext + (info_->sssr() ? 0 : is);
            const shared_ptr<Tensor_<double>> l = l_orthogonal_->data(istate)->at(pos);
            const shared_ptr<Tensor_<double>> t = t_orthogonal_->data(istate)->at(pos);
            const shared_ptr<Tensor_<double>> d = t_orthogonal_->denom(istate)->at(pos);
            for (auto& i3 : t_orthogonal_->virt())
              for (auto& i2 : t_orthogonal_->closed())
                for (auto& i1 : t_orthogonal_->virt())
                  for (auto& i0 : t_orthogonal_->closed()) {
                    if (!t->is_local(i0, i1, i2, i3)) continue;
                    const unique_ptr<double[]> amplitude = t->get_block(i0, i1, i2, i3);
                    const unique_ptr<double[]> lambda    = l->get_block(i0, i1, i2, i3);
                    const unique_ptr<double[]> denom     = d->get_block(i0, i1, i2, i3);
                    unique_ptr<double[]> lcovar          = l->get_block(i0, i1, i2, i3);
                    {
                      const unique_ptr<double[]> lambda2    = l->get_block(i0, i3, i2, i1);
                      sort_indices<0,3,2,1,8,1,-4,1>(lambda2, lcovar, i0.size(), i3.size(), i2.size(), i1.size());
                    }
                    for (size_t j3 = i3.offset(), jall = 0; j3 != i3.offset()+i3.size(); ++j3) {
                      const size_t j3i = j3 - ncore;
                      for (size_t j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) {
                        const size_t j2i = j2 - ncore;
                        for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                          const size_t j1i = j1 - ncore;
                          for (size_t j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++jall) {
                            const size_t j0i = j0 - ncore;
                            const double Lambda = shift2 * lcovar[jall] * amplitude[jall];
                            dshift->element(j0i, j0i) += Lambda;
                            dshift->element(j1i, j1i) -= Lambda;
                            dshift->element(j2i, j2i) += Lambda;
                            dshift->element(j3i, j3i) -= Lambda;
                            nimag[istate] += Lambda;
                            nimag[is] -= Lambda;
                          }
                        }
                      }
                    }
                  }
          }
        }
        timer.tick_print("dshift aibj");
        break;
    }
  }
  dshift->allreduce();
  for (size_t is = 0; is != nstates; ++is)
    for (size_t js = 0; js != nstates; ++js)
      if (!info_->sssr() || is == js) {
        e1->at(is,js)->allreduce();
        e2->at(is,js)->allreduce();
        e3->at(is,js)->allreduce();
        e4->at(is,js)->allreduce();
      }
  mpi__->allreduce(nimag.data(), nimag.size());

  return tie(dshift, e1, e2, e3, e4, nimag);
}

#endif
