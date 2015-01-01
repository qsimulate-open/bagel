//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/smith/moint.h>
#include <src/smith/spinfreebase.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


SpinFreeMethod::SpinFreeMethod(shared_ptr<const SMITH_Info> r) : ref_(r) {
  Timer timer;
  const int max = r->maxtile();
  if (r->ncore() > r->nclosed())
    throw runtime_error("frozen core has been specified but there are not enough closed orbitals");
  closed_ = IndexRange(r->nclosed()-r->ncore(), max, 0, r->ncore());
  active_ = IndexRange(r->nact(),    max, closed_.nblock(),                  r->ncore()+closed_.size());
  virt_   = IndexRange(r->nvirt(),   max, closed_.nblock()+active_.nblock(), r->ncore()+closed_.size()+active_.size());
  all_    = closed_; all_.merge(active_); all_.merge(virt_);
  assert(closed_.size() == r->nclosed() - r->ncore());
  assert(active_.size() == r->nact());
  assert(virt_.size() == r->nvirt());

  rclosed_ = make_shared<const IndexRange>(closed_);
  ractive_ = make_shared<const IndexRange>(active_);
  rvirt_   = make_shared<const IndexRange>(virt_);

  if (ref_->ciwfn()) {
    shared_ptr<const Dvec> dci0 = r->civectors();
    civec_ = dci0->data(ref_->target());
    det_ = civec_->det();

    // length of the ci expansion
    const size_t ci_size = r->civectors()->data(ref_->target())->size();
    ci_ = IndexRange(ci_size, max);
    rci_ = make_shared<const IndexRange>(ci_);
  }

  // f1 tensor.
  {
    vector<IndexRange> o = {all_, all_};
    MOFock fock(ref_, o);
    f1_ = fock.tensor();
    h1_ = fock.h1();
    // canonical orbitals within closed and virtual subspaces
    coeff_ = fock.coeff();
  }

  // for later use
  const int nact = ref_->nact();
  const int nclo = ref_->nclosed();
  auto fockact = make_shared<Matrix>(nact, nact);
  for (auto& i1 : active_)
    for (auto& i0 : active_)
      fockact->copy_block(i0.offset()-nclo, i1.offset()-nclo, i0.size(), i1.size(), f1_->get_block(i0, i1).get());


  // v2 tensor.
  {
    IndexRange occ(closed_);
    occ.merge(active_);
    IndexRange virt(active_);
    virt.merge(virt_);

    vector<IndexRange> o = {occ, virt, occ, virt};
    K2ext v2k(ref_, coeff_, o);
    v2_ = v2k.tensor();
  }

  timer.tick_print("MO integral evaluation");

  // make a ci tensor.
  if (ref_->ciwfn()) {
    vector<IndexRange> o = {ci_};
    // TODO fix later when referece has civec
    Ci dci(ref_, o, civec_);
    rdm0deriv_ = dci.tensor();
  }

  // rdm ci derivatives.
  if (ref_->ciwfn()) {
    shared_ptr<const Dvec> rdm1d = r->rdm1deriv(ref_->target());

    vector<IndexRange> o = {ci_, active_, active_};
    rdm1deriv_ = make_shared<Tensor>(o, false);
    for (auto& i0 : active_) {
      for (auto& i1 : active_) {
        for (auto& ci0 : ci_) {
          const size_t size = i0.size() * i1.size() * ci0.size();
          unique_ptr<double[]> data(new double[size]);
          int iall = 0;
          for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is creation
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
              for (int j2 = ci0.offset(); j2 != ci0.offset()+ci0.size(); ++j2, ++iall)
                // Dvec - first index is annihilation, second is creation (see const_phis_ in fci/determinants.h and knowles_compute.cc)
                data[iall] = rdm1d->data((j1-nclo)+r->nact()*(j0-nclo))->data(j2);
          rdm1deriv_->put_block(data, ci0, i1, i0);
        }
      }
    }
  }

  if (ref_->ciwfn()) {
    shared_ptr<const Dvec> rdm2d = r->rdm2deriv(ref_->target());

    vector<IndexRange> o = {ci_, active_, active_, active_, active_};
    rdm2deriv_ = make_shared<Tensor>(o, false);
    const int nclo = ref_->nclosed();
    for (auto& i0 : active_) {
      for (auto& i1 : active_) {
        for (auto& i2 : active_) {
          for (auto& i3 : active_) {
            for (auto& ci0 : ci_) {
              const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * ci0.size();
              unique_ptr<double[]> data(new double[size]);
              int iall = 0;
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is creation
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
                  for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) // this is creation
                    for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3) // this is annihilation
                      for (int j4 = ci0.offset(); j4 != ci0.offset()+ci0.size(); ++j4, ++iall)
                        data[iall] = rdm2d->data((j3-nclo)+r->nact()*((j2-nclo)+r->nact()*((j1-nclo)+r->nact()*(j0-nclo))))->data(j4);
              rdm2deriv_->put_block(data, ci0, i3, i2, i1, i0);
            }
          }
        }
      }
    }
  }

  if (ref_->ciwfn()) {
    shared_ptr<const Dvec> rdm3d = r->rdm3deriv(ref_->target());
    // RDM4 is contracted a priori by the Fock operator
    shared_ptr<const Dvec> rdm4d = r->rdm4deriv(ref_->target(), fockact);
    assert(rdm3d->ij() == rdm4d->ij());

    vector<IndexRange> o = {ci_, active_, active_, active_, active_, active_, active_};
    rdm3deriv_ = make_shared<Tensor>(o, false);
    rdm4deriv_ = make_shared<Tensor>(o, false);
    const int nclo = ref_->nclosed();
    for (auto& i0 : active_) {
      for (auto& i1 : active_) {
        for (auto& i2 : active_) {
          for (auto& i3 : active_) {
            for (auto& i4 : active_) {
              for (auto& i5 : active_) {
                for (auto& ci0 : ci_) {
                  const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size() * ci0.size();
                  unique_ptr<double[]> data(new double[size]);
                  unique_ptr<double[]> data2(new double[size]);
                  int iall = 0;
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0) // this is  creation
                    for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1) // this is annihilation
                      for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2) // this is creation
                        for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3) // this is annihilation
                          for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4) // this is creation
                            for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5) // this is annhilation
                              for (int j6 = ci0.offset(); j6 != ci0.offset()+ci0.size(); ++j6, ++iall) {
                                data[iall]  = rdm3d->data((j5-nclo)+r->nact()*((j4-nclo)+r->nact()*((j3-nclo)+r->nact()*((j2-nclo)+r->nact()*((j1-nclo)+r->nact()*((j0-nclo)))))))->data(j6);
                                data2[iall] = rdm4d->data((j5-nclo)+r->nact()*((j4-nclo)+r->nact()*((j3-nclo)+r->nact()*((j2-nclo)+r->nact()*((j1-nclo)+r->nact()*((j0-nclo)))))))->data(j6);
                              }
                  rdm3deriv_->put_block(data, ci0, i5, i4, i3, i2, i1, i0);
                  rdm4deriv_->put_block(data2, ci0, i5, i4, i3, i2, i1, i0);
                }
              }
            }
          }
        }
      }
    }
  }

  timer.tick_print("RDM derivative evaluation");

  // rdms.
  if (ref_->ciwfn()) {
    vector<IndexRange> o = {active_, active_};
    rdm1_ = make_shared<Tensor>(o, false);
    const int nclo = ref_->nclosed();
    for (auto& i1 : active_) {
      for (auto& i0 : active_) {
        const size_t size = i0.size() * i1.size();
        unique_ptr<double[]> data(new double[size]);
        int iall = 0;
        for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
          for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
            data[iall] = ref_->rdm1(ref_->target())->element(j0-nclo, j1-nclo);
        rdm1_->put_block(data, i0, i1);
      }
    }
  }
  if (ref_->ciwfn()) {
    vector<IndexRange> o = {active_, active_, active_, active_};
    rdm2_ = make_shared<Tensor>(o, false);
    const int nclo = ref_->nclosed();
    for (auto& i3 : active_) {
      for (auto& i2 : active_) {
        for (auto& i1 : active_) {
          for (auto& i0 : active_) {
            const size_t size = i0.size() * i1.size() * i2.size() * i3.size();
            unique_ptr<double[]> data(new double[size]);
            int iall = 0;
            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                    data[iall] = ref_->rdm2(ref_->target())->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo);
             rdm2_->put_block(data, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  // TODO generic function??
  if (ref_->ciwfn()) {
    {
      vector<IndexRange> o = {active_, active_, active_, active_, active_, active_};
      rdm3_ = make_shared<Tensor>(o, false);
      vector<IndexRange> p = {active_, active_, active_, active_, active_, active_, active_, active_};
      rdm4_ = make_shared<Tensor>(p, false);
    }

    shared_ptr<RDM<3>> rdm3;
    shared_ptr<RDM<4>> rdm4;
    tie(rdm3, rdm4) = ref_->compute_rdm34(ref_->target());

    const int nclo = ref_->nclosed();
    for (auto& i5 : active_)
      for (auto& i4 : active_)
        for (auto& i3 : active_)
          for (auto& i2 : active_)
            for (auto& i1 : active_)
              for (auto& i0 : active_) {
                const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size();
                unique_ptr<double[]> data(new double[size]);
                int iall = 0;
                for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                  for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                    for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                      for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                        for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                          for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                            data[iall] = rdm3->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo);
                rdm3_->put_block(data, i0, i1, i2, i3, i4, i5);
              }
    // TODO there should be a better way of doing this!!!
    for (auto& i7 : active_)
      for (auto& i6 : active_)
        for (auto& i5 : active_)
          for (auto& i4 : active_)
            for (auto& i3 : active_)
              for (auto& i2 : active_)
                for (auto& i1 : active_)
                  for (auto& i0 : active_) {
                    const size_t size = i0.size() * i1.size() * i2.size() * i3.size() * i4.size() * i5.size() * i6.size() * i7.size();
                    unique_ptr<double[]> data(new double[size]);
                    int iall = 0;
                    for (int j7 = i7.offset(); j7 != i7.offset()+i7.size(); ++j7)
                      for (int j6 = i6.offset(); j6 != i6.offset()+i6.size(); ++j6)
                        for (int j5 = i5.offset(); j5 != i5.offset()+i5.size(); ++j5)
                          for (int j4 = i4.offset(); j4 != i4.offset()+i4.size(); ++j4)
                            for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
                              for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
                                for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                                  for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                                    data[iall] = rdm4->element(j0-nclo, j1-nclo, j2-nclo, j3-nclo, j4-nclo, j5-nclo, j6-nclo, j7-nclo);
                rdm4_->put_block(data, i0, i1, i2, i3, i4, i5, i6, i7);
              }


    auto rdm1 = make_shared<RDM<1>>(*ref_->rdm1(ref_->target()));
    auto rdm2 = make_shared<RDM<2>>(*ref_->rdm2(ref_->target()));

    timer.tick_print("RDM evaluation");

    // construct denominator
    denom_ = make_shared<const Denom>(*rdm1, *rdm2, *rdm3, *rdm4, *fockact);

    timer.tick_print("Denominator evaluation");
  }

  // set e0
  e0_ = compute_e0();
}


void SpinFreeMethod::print_iteration() const {
  cout << "      ---- iteration ----" << endl << endl;
  time_ = chrono::high_resolution_clock::now();
}

void SpinFreeMethod::print_iteration(const int i, const double en, const double err) const {
  auto end = chrono::high_resolution_clock::now();
  const double tim = chrono::duration_cast<chrono::milliseconds>(end-time_).count() * 0.001;
  cout << "     " << setw(4) << i << setw(15) << fixed << setprecision(10) << en
                                            << setw(15) << fixed << setprecision(10) << err
                                            << setw(10) << fixed << setprecision(2) << tim << endl;
  time_ = end;
}

void SpinFreeMethod::print_iteration(const bool noconv) const {
  cout << endl << "      -------------------" << endl;
  if (noconv) cout << "      *** Convergence not reached ***" << endl;
  cout << endl;
}

double SpinFreeMethod::compute_e0() const {
  if (ref_->nact() != 0 && !(static_cast<bool>(f1_) && static_cast<bool>(rdm1_)))
    throw logic_error("SpinFreeMethod::compute_e0 was called before f1_ or rdm1_ was computed. Strange.");
  double sum = 0.0;
  for (auto& i1 : active_) {
    for (auto& i0 : active_) {
      const size_t size = i0.size() * i1.size();
      unique_ptr<double[]> fdata = f1_->get_block(i0, i1);
      unique_ptr<double[]> rdata = rdm1_->get_block(i0, i1);
      sum += ddot_(size, fdata, 1, rdata, 1);
    }
  }
  cout << "    - Zeroth order energy: " << setw(20) << setprecision(10) << sum << endl;
  return sum;
}


void SpinFreeMethod::update_amplitude(shared_ptr<Tensor> t, const shared_ptr<Tensor> r, const bool put) {

  // ranks of t and r are assumed to be the same

  // TODO should be parallelized
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          unique_ptr<double[]>       data0 = r->get_block(i0, i1, i2, i3);
          const unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);

          // this is an inverse of the overlap.
          // prefactor of 0.25 included here
          sort_indices<0,3,2,1,2,12,1,12>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                  // note that e0 is cancelled by another term
                  data0[iall] /= (eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
          if (!put) {
            t->add_block(data0, i0, i1, i2, i3);
          } else {
            t->put_block(data0, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i2 : active_) {
    for (auto& i0 : active_) {
      // trans is the transformation matrix
      assert(shalf_xx());
      const int nact = ref_->nact();
      const int nclo = ref_->nclosed();
      unique_ptr<double[]> transp(new double[i0.size()*i2.size()*nact*nact]);
      for (int j2 = i2.offset(), k = 0; j2 != i2.offset()+i2.size(); ++j2)
        for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
          copy_n(shalf_xx()->element_ptr(0,(j0-nclo)+(j2-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

      for (auto& i3 : virt_) {
        for (auto& i1 : virt_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          // data0 is the source area
          unique_ptr<double[]> data0 = r->get_block(i0, i1, i2, i3);
          unique_ptr<double[]> data1(new double[r->get_size(i0, i1, i2, i3)]);
          // sort. Active indices run faster
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i1.size()*i3.size()*nact*nact]);

          // move to orthogonal basis
          dgemm_("N", "N", nact*nact, i1.size()*i3.size(), i0.size()*i2.size(), 1.0, transp, nact*nact, data1, i0.size()*i2.size(),
                                                                                0.0, interm, nact*nact);

          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j02 = 0; j02 != nact*nact; ++j02, ++iall)
                interm[iall] /= e0_ - (denom_xx(j02) + eig_[j3] + eig_[j1]);

          // move back to non-orthogonal basis
          // factor of 0.5 due to the factor in the overlap
          dgemm_("T", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, transp, nact*nact, interm, nact*nact,
                                                                                0.0, data0,  i0.size()*i2.size());

          // sort back to the original order
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
          if (!put) {
            t->add_block(data1, i0, i1, i2, i3);
          } else {
            t->put_block(data1, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i0 : active_) {
    // trans is the transformation matrix
    assert(shalf_x());
    const int nact = ref_->nact();
    const int nclo = ref_->nclosed();
    unique_ptr<double[]> transp(new double[i0.size()*nact]);
    for (int j0 = i0.offset(), k = 0; j0 != i0.offset()+i0.size(); ++j0, ++k)
      copy_n(shalf_x()->element_ptr(0,j0-nclo), nact, transp.get()+nact*k);

    for (auto& i3 : virt_) {
      for (auto& i2 : closed_) {
        for (auto& i1 : virt_) {
          if (!r->get_size_alloc(i2, i3, i0, i1)) continue;
          assert(r->get_size_alloc(i2, i1, i0, i3));
          unique_ptr<double[]>       data0 = r->get_block(i2, i3, i0, i1);
          const unique_ptr<double[]> data1 = r->get_block(i2, i1, i0, i3);
          unique_ptr<double[]> data2(new double[r->get_size(i2, i3, i0, i1)]);
          sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
          sort_indices<2,1,0,3,2,3,1,3>(data1, data2, i2.size(), i1.size(), i0.size(), i3.size());

          // move to orthogonal basis
          unique_ptr<double[]> interm(new double[i1.size()*i2.size()*i3.size()*nact]);
          dgemm_("N", "N", nact, i1.size()*i2.size()*i3.size(), i0.size(), 1.0, transp, nact, data2, i0.size(),
                                                                           0.0, interm, nact);

          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = 0; j0 != nact; ++j0, ++iall)
                  interm[iall] /= e0_ - (denom_x(j0) + eig_[j3] - eig_[j2] + eig_[j1]);

          // move back to non-orthogonal basis
          dgemm_("T", "N", i0.size(), i1.size()*i2.size()*i3.size(), nact, 1.0, transp, nact, interm, nact,
                                                                           0.0, data2,  i0.size());

          if (!put) {
            t->add_block(data2, i0, i1, i2, i3);
          } else {
            t->put_block(data2, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i3 : active_) {
    // trans is the transformation matrix
    assert(shalf_h());
    const int nact = ref_->nact();
    const int nclo = ref_->nclosed();
    unique_ptr<double[]> transp(new double[i3.size()*nact]);
    for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3, ++k)
      copy_n(shalf_h()->element_ptr(0,j3-nclo), nact, transp.get()+nact*k);

    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          if (!r->get_size_alloc(i2, i3, i0, i1)) continue;
          assert(r->get_size_alloc(i0, i3, i2, i1));
          unique_ptr<double[]>       data0 = r->get_block(i2, i3, i0, i1);
          const unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);
          unique_ptr<double[]> data2(new double[r->get_size(i2, i3, i0, i1)]);
          sort_indices<2,3,0,1,0,1,1,1>(data0, data2, i2.size(), i3.size(), i0.size(), i1.size());
          sort_indices<0,3,2,1,2,3,1,3>(data1, data2, i0.size(), i3.size(), i2.size(), i1.size());
          unique_ptr<double[]> interm(new double[i0.size()*i1.size()*i2.size()*nact]);

          // move to orthogonal basis
          dgemm_("N", "T", i0.size()*i1.size()*i2.size(), nact, i3.size(), 1.0, data2, i0.size()*i1.size()*i2.size(), transp, nact,
                                                                           0.0, interm, i0.size()*i1.size()*i2.size());

          size_t iall = 0;
          for (int j3 = 0; j3 != nact; ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                  interm[iall] /= e0_ - (denom_h(j3) - eig_[j2] + eig_[j1] - eig_[j0]);

          // move back to non-orthogonal basis
          dgemm_("N", "N", i0.size()*i1.size()*i2.size(), i3.size(), nact, 1.0, interm, i0.size()*i1.size()*i2.size(), transp, nact,
                                                                           0.0, data2,  i0.size()*i1.size()*i2.size());

          if (!put) {
            t->add_block(data2, i0, i1, i2, i3);
          } else {
            t->put_block(data2, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i3 : active_) {
    for (auto& i1 : active_) {
      assert(shalf_hh());
      const int nact = ref_->nact();
      const int nclo = ref_->nclosed();
      unique_ptr<double[]> transp(new double[i1.size()*i3.size()*nact*nact]);
      for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
        for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++k)
          copy_n(shalf_hh()->element_ptr(0,(j1-nclo)+(j3-nclo)*nact), nact*nact, transp.get()+nact*nact*k);

      for (auto& i2 : closed_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          // data0 is the source area
          unique_ptr<double[]> data0 = r->get_block(i0, i1, i2, i3);
          unique_ptr<double[]> data1(new double[r->get_size(i0, i1, i2, i3)]);
          // sort. Active indices run slower
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i1.size(), i2.size(), i3.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i0.size()*i2.size()*nact*nact]);

          // move to orthogonal basis
          dgemm_("N", "T", i0.size()*i2.size(), nact*nact, i1.size()*i3.size(), 1.0, data1, i0.size()*i2.size(), transp, nact*nact,
                                                                                0.0, interm, i0.size()*i2.size());

          size_t iall = 0;
          for (int j13 = 0; j13 != nact*nact; ++j13)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                interm[iall] /= e0_ - (denom_hh(j13) - eig_[j2] - eig_[j0]);

          // move back to non-orthogonal basis
          // factor of 0.5 due to the factor in the overlap
          dgemm_("N", "N", i0.size()*i2.size(), i1.size()*i3.size(), nact*nact, 0.5, interm, i0.size()*i2.size(), transp, nact*nact,
                                                                                0.0, data0,  i0.size()*i2.size());

          // sort back to the original order
          sort_indices<0,2,1,3,0,1,1,1>(data0, data1, i0.size(), i2.size(), i1.size(), i3.size());
          if (!put) {
            t->add_block(data1, i0, i1, i2, i3);
          } else {
            t->put_block(data1, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i3 : active_) {
    for (auto& i2 : active_) {
      assert(shalf_xh());
      const int nact = ref_->nact();
      const int nclo = ref_->nclosed();
      unique_ptr<double[]> transp(new double[i2.size()*i3.size()*nact*nact*4]);
      for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
        for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++k) {
          copy_n(shalf_xh()->element_ptr(0,             (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*k);
          copy_n(shalf_xh()->element_ptr(0, nact*nact + (j2-nclo)+(j3-nclo)*nact), nact*nact*2, transp.get()+nact*nact*2*(k+i2.size()*i3.size()));
        }

      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          const size_t blocksize = r->get_size_alloc(i2, i3, i0, i1);
          if (!blocksize) continue;
          assert(blocksize == r->get_size_alloc(i0, i3, i2, i1));
          unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
          unique_ptr<double[]> data1 = r->get_block(i0, i3, i2, i1);

          unique_ptr<double[]> data2(new double[blocksize*2]);
          // sort. Active indices run slower
          sort_indices<2,3,0,1,0,1,1,1>(data0.get(), data2.get()          , i2.size(), i3.size(), i0.size(), i1.size());
          sort_indices<0,3,2,1,0,1,1,1>(data1.get(), data2.get()+blocksize, i0.size(), i3.size(), i2.size(), i1.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i0.size()*i1.size()*nact*nact*2]);

          // move to orthogonal basis
          dgemm_("N", "T", i0.size()*i1.size(), nact*nact*2, i2.size()*i3.size()*2, 1.0, data2,  i0.size()*i1.size(), transp, nact*nact*2,
                                                                                    0.0, interm, i0.size()*i1.size());

          size_t iall = 0;
          for (int j23 = 0; j23 != nact*nact*2; ++j23)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
              for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                interm[iall] /= e0_ - (denom_xh(j23) + eig_[j1] - eig_[j0]);

          // move back to non-orthogonal basis
          dgemm_("N", "N", i0.size()*i1.size(), i2.size()*i3.size()*2, nact*nact*2, 1.0, interm, i0.size()*i1.size(), transp, nact*nact*2,
                                                                                    0.0, data2,  i0.size()*i1.size());

          // sort back to the original order
          copy_n(data2.get(), blocksize, data0.get());
          sort_indices<2,1,0,3,0,1,1,1>(data2.get()+blocksize, data1.get(), i0.size(), i1.size(), i2.size(), i3.size());
          if (!put) {
            t->add_block(data0, i0, i1, i2, i3);
            t->add_block(data1, i2, i1, i0, i3);
          } else {
            t->put_block(data0, i0, i1, i2, i3);
            t->put_block(data1, i2, i1, i0, i3);
          }
        }
      }
    }
  }
  for (auto& i3 : active_) {
    for (auto& i2 : active_) {
      for (auto& i0 : active_) {
        assert(shalf_xhh());
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        unique_ptr<double[]> transp(new double[i0.size()*i2.size()*i3.size()*nact*nact*nact]);
        for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
          for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
              copy_n(shalf_xhh()->element_ptr(0,j0-nclo+nact*(j2-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

        for (auto& i1 : virt_) {
          // if this block is not included in the current wave function, skip it
          const size_t blocksize = r->get_size_alloc(i2, i3, i0, i1);
          if (!blocksize) continue;
          // data0 is the source area
          unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
          unique_ptr<double[]> data1(new double[blocksize]);
          // sort. Active indices run slower
          sort_indices<3,2,0,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i1.size()*nact*nact*nact]);

          // move to orthogonal basis
          dgemm_("N", "T", i1.size(), nact*nact*nact, i0.size()*i2.size()*i3.size(), 1.0, data1,  i1.size(), transp, nact*nact*nact,
                                                                                     0.0, interm, i1.size());

          size_t iall = 0;
          for (int j123 = 0; j123 != nact*nact*nact; ++j123)
            for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1, ++iall)
              interm[iall] /= e0_ - (denom_xhh(j123) + eig_[j1]);

          // move back to non-orthogonal basis
          dgemm_("N", "N", i1.size(), i0.size()*i2.size()*i3.size(), nact*nact*nact, 1.0, interm, i1.size(), transp, nact*nact*nact,
                                                                                     0.0, data0,  i1.size());

          // sort back to the original order
          sort_indices<1,0,2,3,0,1,1,1>(data0, data1, i1.size(), i0.size(), i2.size(), i3.size());
          if (!put) {
            t->add_block(data1, i0, i1, i2, i3);
          } else {
            t->put_block(data1, i0, i1, i2, i3);
          }
        }
      }
    }
  }
  for (auto& i3 : active_) {
    for (auto& i1 : active_) {
      for (auto& i0 : active_) {
        assert(shalf_xxh());
        const int nact = ref_->nact();
        const int nclo = ref_->nclosed();
        unique_ptr<double[]> transp(new double[i0.size()*i1.size()*i3.size()*nact*nact*nact]);
        for (int j3 = i3.offset(), k = 0; j3 != i3.offset()+i3.size(); ++j3)
          for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
            for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++k)
              copy_n(shalf_xxh()->element_ptr(0,j0-nclo+nact*(j1-nclo+nact*(j3-nclo))), nact*nact*nact, transp.get()+nact*nact*nact*k);

        for (auto& i2 : closed_) {
          // if this block is not included in the current wave function, skip it
          const size_t blocksize = r->get_size_alloc(i2, i3, i0, i1);
          if (!blocksize) continue;
          // data0 is the source area
          unique_ptr<double[]> data0 = r->get_block(i2, i3, i0, i1);
          unique_ptr<double[]> data1(new double[blocksize]);
          // sort. Active indices run slower
          sort_indices<0,2,3,1,0,1,1,1>(data0, data1, i2.size(), i3.size(), i0.size(), i1.size());
          // intermediate area
          unique_ptr<double[]> interm(new double[i2.size()*nact*nact*nact]);

          // move to orthogonal basis
          dgemm_("N", "T", i2.size(), nact*nact*nact, i0.size()*i1.size()*i3.size(), 1.0, data1,  i2.size(), transp, nact*nact*nact,
                                                                                     0.0, interm, i2.size());

          size_t iall = 0;
          for (int j013 = 0; j013 != nact*nact*nact; ++j013)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2, ++iall)
              interm[iall] /= e0_ - (denom_xxh(j013) - eig_[j2]);

          // move back to non-orthogonal basis
          dgemm_("N", "N", i2.size(), i0.size()*i1.size()*i3.size(), nact*nact*nact, 1.0, interm, i2.size(), transp, nact*nact*nact,
                                                                                     0.0, data0,  i2.size());

          // sort back to the original order
          sort_indices<1,2,0,3,0,1,1,1>(data0, data1, i2.size(), i0.size(), i1.size(), i3.size());
          if (!put) {
            t->add_block(data1, i0, i1, i2, i3);
          } else {
            t->put_block(data1, i0, i1, i2, i3);
          }
        }
      }
    }
  }
}
