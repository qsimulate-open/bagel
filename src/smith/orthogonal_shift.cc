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
    switch(iext) {
      case Excitations::arbs:
        for (auto& i3 : virt_)
          for (auto& i1 : virt_)
            for (auto& i0o : interm_[iext]) {
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
              for (auto& i0o : interm_[iext]) {
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
              for (auto& i0o : interm_[iext]) {
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
            for (auto& i0o : interm_[iext]) {
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
            for (auto& i0o : interm_[iext]) {
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
          for (auto& i0o : interm_[iext]) {
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
          for (auto& i0o : interm_[iext]) {
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
    

tuple<shared_ptr<Matrix>,shared_ptr<Matrix>> Orthogonal_Basis::get_zX(const int iext, shared_ptr<const Matrix> largey, shared_ptr<const Matrix> smallzin) const {
  shared_ptr<Matrix> smallz = smallzin->copy();
  shared_ptr<Matrix> largex = largey->clone();

  const size_t interm_size = shalf_[iext]->ndim();
  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
    for (size_t j1o = 0; j1o != interm_size; ++j1o) {
      if (j1o == j0o) continue;
      const double fdiff = (phi_[iext])[j1o] - (phi_[iext])[j0o];
      smallz->element(j1o, j0o) = fabs(fdiff) > 1.0e-12 ? -0.5 * (largey->element(j1o, j0o) - largey->element(j0o, j1o)) / fdiff : 0.0;
    }
  }
  for (size_t j0o = 0; j0o != interm_size; ++j0o) {
    for (size_t j1o = 0; j1o != interm_size; ++j1o) {
      largex->element(j1o, j0o) = 0.25 * (largey->element(j1o, j0o) + 2.0 * smallz->element(j1o, j0o) * (phi_[iext])[j1o])
                                + 0.25 * (largey->element(j0o, j1o) + 2.0 * smallz->element(j0o, j1o) * (phi_[iext])[j0o]);
    }
  }
  
  return tie(smallz, largex);
}


tuple<shared_ptr<Matrix>,shared_ptr<Vec<double>>,shared_ptr<VecRDM<1>>,shared_ptr<VecRDM<2>>,shared_ptr<VecRDM<3>>,shared_ptr<VecRDM<3>>,vector<double>>
Orthogonal_Basis::make_d2_imag(shared_ptr<const Orthogonal_Basis> lambda) const {
  auto dshift = make_shared<Matrix>(norb_-ncore_, norb_-ncore_);
  auto e0 = make_shared<Vec<double>>();
  auto e1 = make_shared<VecRDM<1>>();
  auto e2 = make_shared<VecRDM<2>>();
  auto e3 = make_shared<VecRDM<3>>();
  auto e4 = make_shared<VecRDM<3>>();
  vector<double> nimag;
  nimag.resize(nstates_);

  for (size_t is = 0; is != nstates_; ++is)
    for (size_t js = 0; js != nstates_; ++js)
      if (!sssr_ || is == js) {
        auto e0temp = make_shared<double>(0.0);
        auto e1temp = make_shared<RDM<1>>(nact_);
        auto e2temp = make_shared<RDM<2>>(nact_);
        auto e3temp = make_shared<RDM<3>>(nact_);
        auto e4temp = make_shared<RDM<3>>(nact_);

        e0->emplace(is, js, e0temp);
        e1->emplace(is, js, e1temp);
        e2->emplace(is, js, e2temp);
        e3->emplace(is, js, e3temp);
        e4->emplace(is, js, e4temp);
      }

//  const double shift2 = shift_ * shift_;

  for (int istate = 0; istate != nstates_; ++istate) {
    shared_ptr<MultiTensor_<double>> tcovar = get_contravariant(istate);
    shared_ptr<MultiTensor_<double>> lcovar = lambda->get_contravariant(istate);
    for (int iext = Excitations::arbs; iext != Excitations::total; ++iext) {
      const shared_ptr<Tensor_<double>> ttensor = data_[istate]->at(iext);
      const shared_ptr<Tensor_<double>> dtensor = denom_[istate]->at(iext);
      const shared_ptr<Tensor_<double>> ltensor = lambda->data(istate)->at(iext);
      const size_t interm_size = (iext == Excitations::aibj) ? 0 : shalf_[iext]->ndim();
      auto largey = make_shared<Matrix>(interm_size, interm_size);
      auto largeq = make_shared<Matrix>(interm_size, interm_size);
      auto smallz = make_shared<Matrix>(interm_size, interm_size);
      shared_ptr<Matrix> largex;
      switch(iext) {
        // arbs is already wrong. should debug.............. 
        case Excitations::arbs:
#if 0
          for (auto& i3 : virt_)
            for (auto& i1 : virt_)
              for (auto& i0o : interm_[iext]) {
                if (!ttensor->is_local(i0o, i1, i3)) continue;
                const unique_ptr<double[]> t0data = ttensor->get_block(i0o, i1, i3);
                const unique_ptr<double[]> d0data = dtensor->get_block(i0o, i1, i3);
                const unique_ptr<double[]> l0data = ltensor->get_block(i0o, i1, i3);
                for (size_t j3 = i3.offset(), i0all = 0; j3 != i3.offset()+i3.size(); ++j3) {
                  for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                    for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++i0all) {
                      const double Lambda = shift2 * l0data[i0all] * t0data[i0all];
                      dshift->element(j1-ncore_, j1-ncore_) -= Lambda;
                      dshift->element(j3-ncore_, j3-ncore_) -= Lambda;
                      smallz->element(j0o, j0o) -= Lambda;
                      nimag[istate] -= Lambda;
                    }
                  }
                }
                for (auto& i1o : interm_[iext]) {
                  const unique_ptr<double[]> t1data = ttensor->get_block(i1o, i1, i3);
                  const unique_ptr<double[]> d1data = dtensor->get_block(i1o, i1, i3);
                  const unique_ptr<double[]> l1data = ltensor->get_block(i1o, i1, i3);
                  for (size_t j3 = i3.offset(), i0all = 0; j3 != i3.offset()+i3.size(); ++j3) {
                    for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                      for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++i0all) {
                        for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                          const size_t i1all = j1o-i1o.offset() + i1o.size() * (j1-i1.offset() + i1.size() * (j3-i3.offset()));
                          const double lt = l0data[i0all] * t1data[i1all] * shift2 * d0data[i0all];
                          const double tl = l1data[i1all] * t0data[i0all] * shift2 * d1data[i1all];
                          largey->element(j0o, j1o) += l0data[i0all] * t1data[i1all] * shift2 * (d0data[i0all] - d1data[i1all]);
                          largeq->element(j0o, j1o) += (lt + tl);
                        }
                      }
                    }
                  }
                }
              }
          tie(smallz, largex) = get_zX(iext, largey, smallz);
          // will tensor parallelize it later...
          for (size_t is = 0; is != nstates_; ++is) {
            for (size_t js = 0; js != nstates_; ++js) {
              if (sssr_ && is != js) continue;
              shared_ptr<RDM<1>> rdm1;
              shared_ptr<RDM<2>> rdm2;
              shared_ptr<RDM<3>> rdm3;
              shared_ptr<RDM<4>> rdm4;
              tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(is, js);
              {
                shared_ptr<Matrix> Vmat = shalf_[iext]->get_submatrix(0, is * nact_ * nact_, interm_size, nact_ * nact_);
                auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
                auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));
                auto Rmat = make_shared<RDM<2>>(nact_);
                // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
                for (size_t j0 = 0; j0 != nact_; ++j0) {
                  for (size_t j1 = 0; j1 != nact_; ++j1) {
                    for (size_t j2 = 0; j2 != nact_; ++j2) {
                      for (size_t j3 = 0; j3 != nact_; ++j3) {
                        for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                          const double VrsO = shalf_[iext]->element(j0o, j2 + j3 * nact_ + js * nact_ * nact_);
                          Rmat->element(j0, j1, j2, j3) += Qmat->element(j0o, j0 + j1 * nact_) * VrsO;
                          e2->at(is, js)->element(j0, j2, j1, j3) += Pmat->element(j0o, j0 + j1 * nact_) * VrsO;
                        }
                      }
                    }
                  }
                }
                // (3) form d^{(2)}_{rs} and e3
                for (size_t j0 = 0; j0 != nact_; ++j0)
                  for (size_t j1 = 0; j1 != nact_; ++j1)
                    for (size_t j2 = 0; j2 != nact_; ++j2)
                      for (size_t j3 = 0; j3 != nact_; ++j3)
                        for (size_t j4 = 0; j4 != nact_; ++j4) {
                          const size_t j4i = j4 + nclo_;
                          for (size_t j5 = 0; j5 != nact_; ++j5) {
                            const size_t j5i = j5 + nclo_;
                            dshift->element(j4i, j5i) += Rmat->element(j0, j1, j2, j3) * rdm3->element(j0, j2, j1, j3, j4, j5);
                            e3->at(is, js)->element(j0, j2, j1, j3, j4, j5) += 2.0 * Rmat->element(j0, j1, j2, j3) * fockact_->element(j4, j5);
                          }
                        }
              }
            }
          }
#endif
          break;
        case Excitations::arbi:
#if 0
          for (auto& i3 : virt_)
            for (auto& i2 : closed_)
              for (auto& i1 : virt_)
                for (auto& i0o : interm_[iext]) {
                  if (!ttensor->is_local(i0o, i1, i2, i3)) continue;
//                  const unique_ptr<double[]> t0data = ttensor->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> t0data = ttensor->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> d0data = dtensor->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> l0data = ltensor->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> lcdata = lcovar->at(iext)->get_block(i0o, i1, i2, i3);
                  const unique_ptr<double[]> tcdata = tcovar->at(iext)->get_block(i0o, i1, i2, i3);
                  for (size_t j3 = i3.offset(), i0all = 0; j3 != i3.offset()+i3.size(); ++j3) {
                    for (size_t j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) {
                      for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                        for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++i0all) {
                          const double Lambda = shift2 * lcdata[i0all] * t0data[i0all];
                          dshift->element(j1-ncore_, j1-ncore_) -= Lambda;
                          dshift->element(j3-ncore_, j3-ncore_) -= Lambda;
                          dshift->element(j2-ncore_, j2-ncore_) += Lambda;
                          smallz->element(j0o, j0o) -= Lambda;
                          nimag[istate] -= Lambda;
                        }
                      }
                    }
                  }
                  for (auto& i1o : interm_[iext]) {
//                    const unique_ptr<double[]> t1data = ttensor->get_block(i1o, i1, i2, i3);
                    const unique_ptr<double[]> t1data = ttensor->get_block(i1o, i1, i2, i3);
                    const unique_ptr<double[]> d1data = dtensor->get_block(i1o, i1, i2, i3);
                    const unique_ptr<double[]> l1data = ltensor->get_block(i1o, i1, i2, i3);
                    for (size_t j3 = i3.offset(), i0all = 0; j3 != i3.offset()+i3.size(); ++j3) {
                      for (size_t j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) {
                        for (size_t j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) {
                          for (size_t j0o = i0o.offset(); j0o != i0o.offset()+i0o.size(); ++j0o, ++i0all) {
                            for (size_t j1o = i1o.offset(); j1o != i1o.offset()+i1o.size(); ++j1o) {
                              const size_t i1all = j1o-i1o.offset() + i1o.size() * (j1-i1.offset() + i1.size() * (j2 - i2.offset() + i2.size() * (j3-i3.offset())));
                              const double lt = lcdata[i0all] * t1data[i1all] * shift2 * d0data[i0all];
                              const double tl = l1data[i1all] * tcdata[i0all] * shift2 * d1data[i1all];
                              largey->element(j0o, j1o) += lcdata[i0all] * t1data[i1all] * shift2 * (d0data[i0all] - d1data[i1all]);
                              largeq->element(j0o, j1o) += (lt + tl);
                            }
                          }
                        }
                      }
                    }
                  }
                }
          tie(smallz, largex) = get_zX(iext, largey, smallz);
          for (size_t is = 0; is != nstates_; ++is) {
            for (size_t js = 0; js != nstates_; ++js) {
              if (sssr_ && is != js) continue;
              shared_ptr<RDM<1>> rdm1;
              shared_ptr<RDM<2>> rdm2;
              shared_ptr<RDM<3>> rdm3;
              shared_ptr<RDM<4>> rdm4;
              tie(rdm1, rdm2, rdm3, rdm4) = feed_rdm(js, is);
              {
                shared_ptr<Matrix> Vmat = shalf_[iext]->get_submatrix(0, is * nact_, interm_size, nact_);
                auto Qmat = make_shared<Matrix>((*smallz) % (*Vmat));
                auto Pmat = make_shared<Matrix>((-2.0 * (*largex) + (*largeq)) % (*Vmat));
                auto Rmat = make_shared<RDM<1>>(nact_);
                // (2) form R_{tv,uw} = Q^U_{tv} V^{U}_{uw}
                for (size_t j0 = 0; j0 != nact_; ++j0) {
                  for (size_t j1 = 0; j1 != nact_; ++j1) {
                    for (size_t j0o = 0; j0o != interm_size; ++j0o) {
                      const double VrO = shalf_[iext]->element(j0o, j1 + js * nact_);
                      Rmat->element(j0, j1) += Qmat->element(j0o, j0) * VrO;
                      e1->at(is, js)->element(j0, j1) += Pmat->element(j0o, j0) * VrO;
                    }
                  }
                }
                // (3) form d^{(2)}_{rs} and e3
                for (size_t j0 = 0; j0 != nact_; ++j0)
                  for (size_t j1 = 0; j1 != nact_; ++j1)
                    for (size_t j2 = 0; j2 != nact_; ++j2) {
                      const size_t j2i = j2 + nclo_;
                      for (size_t j3 = 0; j3 != nact_; ++j3) {
                        const size_t j3i = j3 + nclo_;
                        dshift->element(j2i, j3i) += Rmat->element(j0, j1) * rdm2->element(j0, j1, j2, j3);
                        e2->at(is, js)->element(j0, j1, j2, j3) += 2.0 * Rmat->element(j0, j1) * fockact_->element(j2, j3);
                      }
                    }
              }
            }
          }
#endif
          break;
        case Excitations::airj:
          break;
        case Excitations::risj:
          break;
        case Excitations::airs:
          break;
        case Excitations::arst:
          break;
        case Excitations::rist:
          break;
        case Excitations::aibj:
          break;
      }
    }

  }

  return tie(dshift, e0, e1, e2, e3, e4, nimag);
}

#endif
