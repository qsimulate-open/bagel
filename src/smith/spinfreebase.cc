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
#include <src/smith/smith_util.h>

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
    core_energy_ = fock.core_energy();
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
    IndexRange occ(closed_);  occ.merge(active_);
    IndexRange virt(active_); virt.merge(virt_);

    // in the case of MRCI, we need to include all sectors
    if (to_upper(ref_->method()) == "MRCI") {
      occ = all_;
      virt = all_;
    }

    vector<IndexRange> o = {occ, virt, occ, virt};
    K2ext v2k(ref_, coeff_, o);
    v2_ = v2k.tensor();
  }

  timer.tick_print("MO integral evaluation");

  // rdm ci derivatives. Only for gradient computations
  if (ref_->ciwfn() && ref_->grad()) {
    shared_ptr<Dvec> rdm0d = make_shared<Dvec>(civec_, 1);
    shared_ptr<Dvec> rdm1d = r->rdm1deriv(ref_->target());
    shared_ptr<Dvec> rdm2d = r->rdm2deriv(ref_->target());
    shared_ptr<Dvec> rdm3d, rdm4d;
    // RDM4 is contracted a priori by the Fock operator
    tie(rdm3d, rdm4d) = r->rdm34deriv(ref_->target(), fockact);
    assert(rdm3d->ij() == rdm4d->ij());

    vector<IndexRange> o1 = {ci_};
    vector<IndexRange> o3 = {ci_, active_, active_};
    vector<IndexRange> o5 = {ci_, active_, active_, active_, active_};
    vector<IndexRange> o7 = {ci_, active_, active_, active_, active_, active_, active_};
    rdm0deriv_ = make_shared<Tensor>(o1);
    rdm1deriv_ = make_shared<Tensor>(o3);
    rdm2deriv_ = make_shared<Tensor>(o5);
    rdm3deriv_ = make_shared<Tensor>(o7);
    rdm4deriv_ = make_shared<Tensor>(o7);

    const int nclo = ref_->nclosed();
    vector<int> inpoff1(1,0);
    vector<int> inpoff3(2,nclo); inpoff3.push_back(0);
    vector<int> inpoff5(4,nclo); inpoff5.push_back(0);
    vector<int> inpoff7(6,nclo); inpoff7.push_back(0);

    const btas::CRange<1> range1(rdm1d->extent(0)*rdm1d->extent(1));
    const btas::CRange<3> range3(rdm1d->extent(0)*rdm1d->extent(1), nact, nact);
    const btas::CRange<5> range5(rdm2d->extent(0)*rdm2d->extent(1), nact, nact, nact, nact);
    const btas::CRange<7> range7(rdm3d->extent(0)*rdm3d->extent(1), nact, nact, nact, nact, nact, nact);

    rdm0d->resize(range1);
    rdm1d->resize(range3);
    rdm2d->resize(range5);
    rdm3d->resize(range7);
    rdm4d->resize(range7);
    fill_block<1>(rdm0deriv_, rdm0d, inpoff1, o1);
    fill_block<3>(rdm1deriv_, rdm1d, inpoff3, o3);
    fill_block<5>(rdm2deriv_, rdm2d, inpoff5, o5);
    fill_block<7>(rdm3deriv_, rdm3d, inpoff7, o7);
    fill_block<7>(rdm4deriv_, rdm4d, inpoff7, o7);

    timer.tick_print("RDM derivative evaluation");
  }


  // rdms.
  if (ref_->ciwfn()) {
    const int nclo = ref_->nclosed();
    const int nstates = ref_->ciwfn()->nstates();
    rdm1all_ = make_shared<Vec<Tensor>>();
    rdm2all_ = make_shared<Vec<Tensor>>();
    rdm3all_ = make_shared<Vec<Tensor>>();
    rdm4all_ = make_shared<Vec<Tensor>>();

    // TODO this can be reduced to half by bra-ket symmetry
    for (int ist = 0; ist != nstates; ++ist) {
      for (int jst = 0; jst != nstates; ++jst) {

        auto rdm1t = make_shared<Tensor>(vector<IndexRange>(2,active_));
        auto rdm2t = make_shared<Tensor>(vector<IndexRange>(4,active_));
        auto rdm3t = make_shared<Tensor>(vector<IndexRange>(6,active_));
        auto rdm4t = make_shared<Tensor>(vector<IndexRange>(8,active_));

        shared_ptr<const RDM<1>> rdm1;
        shared_ptr<const RDM<2>> rdm2;
        shared_ptr<const RDM<3>> rdm3;
        shared_ptr<const RDM<4>> rdm4;
        tie(rdm1, rdm2) = ref_->rdm12(jst, ist);
        tie(rdm3, rdm4) = ref_->rdm34(jst, ist);

        fill_block<2>(rdm1t, rdm1, vector<int>(2,nclo), vector<IndexRange>(2,active_));
        fill_block<4>(rdm2t, rdm2, vector<int>(4,nclo), vector<IndexRange>(4,active_));
        fill_block<6>(rdm3t, rdm3, vector<int>(6,nclo), vector<IndexRange>(6,active_));
        fill_block<8>(rdm4t, rdm4, vector<int>(8,nclo), vector<IndexRange>(8,active_));

        rdm1all_->emplace(jst, ist, rdm1t);
        rdm2all_->emplace(jst, ist, rdm2t);
        rdm3all_->emplace(jst, ist, rdm3t);
        rdm4all_->emplace(jst, ist, rdm4t);

        if (ist == jst)
          // construct denominator
          denom_.push_back(make_shared<Denom>(*rdm1, *rdm2, *rdm3, *rdm4, *fockact));
      }
    }
    timer.tick_print("RDM + denominator evaluation");
  }

  // set e0
  e0_ = compute_e0();
}


void SpinFreeMethod::set_rdm(const int ist, const int jst) {
  rdm1_ = rdm1all_->at(ist, jst);
  rdm2_ = rdm2all_->at(ist, jst);
  rdm3_ = rdm3all_->at(ist, jst);
  rdm4_ = rdm4all_->at(ist, jst);
}


void SpinFreeMethod::print_iteration() {
  cout << "      ---- iteration ----" << endl << endl;
}


void SpinFreeMethod::print_iteration(const int i, const double en, const double err, const double tim, const int ist) {
  cout << "     " << setw(4) << i;
  if (ist >= 0)
    cout << setw(4) << ist;
  cout << setw(15) << fixed << setprecision(8) << en << setw(15) << fixed << setprecision(8) << err
                                                     << setw(10) << fixed << setprecision(2) << tim << endl;
}


void SpinFreeMethod::print_iteration(const bool noconv) {
  cout << endl << "      -------------------" << endl;
  if (noconv) cout << "      *** Convergence not reached ***" << endl;
  cout << endl;
}


double SpinFreeMethod::compute_e0() {
  assert(!!f1_);
  const size_t nstates = ref_->ciwfn()->nstates();
  double sum = 0.0;
  for (int ist = 0; ist != nstates; ++ist) {
    set_rdm(ist, ist);
    assert(!!rdm1_);
    for (auto& i1 : active_) {
      for (auto& i0 : active_) {
        const size_t size = i0.size() * i1.size();
        unique_ptr<double[]> fdata = f1_->get_block(i0, i1);
        unique_ptr<double[]> rdata = rdm1_->get_block(i0, i1);
        sum += ddot_(size, fdata, 1, rdata, 1);
      }
    }
  }
  sum /= nstates;
  cout << "    - Zeroth order energy: " << setw(20) << setprecision(10) << sum << endl;
  return sum;
}


// local function to compress the following
void SpinFreeMethod::loop_over(std::function<void(const Index&, const Index&, const Index&, const Index&)> func) const {
  for (auto& i3 : virt_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i2 : active_)
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i0 : active_)
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_) {
          func(i0, i1, i2, i3);
          func(i2, i1, i0, i3);
        }
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i0 : active_)
        for (auto& i2 : closed_)
          func(i0, i1, i2, i3);
}

shared_ptr<Tensor> SpinFreeMethod::init_amplitude() const {
  shared_ptr<Tensor> out = v2_->clone();
  auto put = [this, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = v2_->get_size_alloc(i0, i1, i2, i3);
    unique_ptr<double[]> buf(new double[size]);
    fill_n(buf.get(), size, 0.0);
    out->put_block(buf, i0, i1, i2, i3);
  };
  loop_over(put);
  return out;
}


shared_ptr<Tensor> SpinFreeMethod::init_residual() const {
  shared_ptr<Tensor> out = v2_->clone();
  auto put = [this, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = v2_->get_size_alloc(i2, i3, i0, i1);
    unique_ptr<double[]> buf(new double[size]);
    fill_n(buf.get(), size, 0.0);
    out->put_block(buf, i2, i3, i0, i1);
  };
  loop_over(put);
  return out;
}


double SpinFreeMethod::dot_product_transpose(shared_ptr<const MultiTensor> r, shared_ptr<const MultiTensor> t2) const {
  assert(r->nref() == t2->nref());
  double out = 0.0;
  for (int i = 0; i != r->nref(); ++i)
    out += r->fac(i) * t2->fac(i)
         + dot_product_transpose(r->at(i), t2->at(i));
  return out;
}


double SpinFreeMethod::dot_product_transpose(shared_ptr<const Tensor> r, shared_ptr<const Tensor> t2) const {
  double out = 0.0;
  auto prod = [this, &r, &t2, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = r->get_size_alloc(i2, i3, i0, i1);
    if (size != 0) {
      unique_ptr<double[]> tmp0 = t2->get_block(i0, i1, i2, i3);
      unique_ptr<double[]> tmp1(new double[size]);
      sort_indices<2,3,0,1,0,1,1,1>(tmp0.get(), tmp1.get(), i0.size(), i1.size(), i2.size(), i3.size());

      out += blas::dot_product(tmp1.get(), size, r->get_block(i2, i3, i0, i1).get());
    }
  };
  loop_over(prod);
  return out;
}


void SpinFreeMethod::diagonal(shared_ptr<Tensor> r, shared_ptr<const Tensor> t) const {
  assert(to_upper(ref_->method()) == "CASPT2");
  for (auto& i3 : virt_) {
    for (auto& i2 : closed_) {
      for (auto& i1 : virt_) {
        for (auto& i0 : closed_) {
          // if this block is not included in the current wave function, skip it
          if (!r->get_size_alloc(i0, i1, i2, i3)) continue;
          unique_ptr<double[]>       data0 = t->get_block(i0, i1, i2, i3);
          const unique_ptr<double[]> data1 = t->get_block(i0, i3, i2, i1);

          sort_indices<0,3,2,1,8,1,-4,1>(data1, data0, i0.size(), i3.size(), i2.size(), i1.size());
          size_t iall = 0;
          for (int j3 = i3.offset(); j3 != i3.offset()+i3.size(); ++j3)
            for (int j2 = i2.offset(); j2 != i2.offset()+i2.size(); ++j2)
              for (int j1 = i1.offset(); j1 != i1.offset()+i1.size(); ++j1)
                for (int j0 = i0.offset(); j0 != i0.offset()+i0.size(); ++j0, ++iall)
                  // note that e0 is cancelled by another term
                  data0[iall] *= -(eig_[j0] + eig_[j2] - eig_[j3] - eig_[j1]);
          r->add_block(data0, i0, i1, i2, i3);
        }
      }
    }
  }
}
