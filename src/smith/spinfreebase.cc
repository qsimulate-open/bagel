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


template<int N>
static void fill_block(shared_ptr<Tensor> target, shared_ptr<const btas::TensorN<double,N>> input,
                       const vector<int> inpoffsets, const vector<IndexRange> ranges) {
  assert(input->range().ordinal().contiguous());
  assert(target->rank() == input->range().rank() && target->rank() > 0);
  const int rank = target->rank();

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  LoopGenerator gen(ranges);
  vector<vector<Index>> loop = gen.block_loop();
  for (auto& indices : loop) {
    assert(indices.size() == rank);

    const size_t buffersize = accumulate(indices.begin(), indices.end(), 1ul, prod);
    unique_ptr<double[]> buffer(new double[buffersize]);
    vector<size_t> stride;
    for (auto i = indices.begin(); i != indices.end(); ++i) {
      auto ii = i; ++ii;
      stride.push_back(accumulate(ii, indices.end(), 1ul, prod));
    }

    vector<size_t> extent(rank);
    auto e = extent.rbegin();
    for (int i = 0; i != rank; ++i)
      *e++ = input->extent(i); 

    vector<size_t> stride_target;
    for (auto i = extent.begin(); i != extent.end(); ++i) {
      auto ii = i; ++ii;
      stride_target.push_back(accumulate(ii, extent.end(), 1ul, multiplies<size_t>()));
    }

    const size_t backsize = indices.back().size();
    for (size_t n = 0; n != buffersize; n += backsize) {
      size_t offset = 0lu;
      size_t tmp = n;
      for (int i = 0; i != rank; ++i) {
        offset += (tmp / stride[i] + indices[i].offset() - inpoffsets[i]) * stride_target[i]; 
        tmp = n % stride[i];
      }
      copy_n(input->data()+offset, backsize, buffer.get()+n);
    }

    target->put_block(buffer, vector<Index>(indices.rbegin(), indices.rend()));
  }
}


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
    vector<IndexRange> o = {ci_};
    Ci dci(ref_, o, civec_);
    rdm0deriv_ = dci.tensor();

    shared_ptr<Dvec> rdm1d = r->rdm1deriv(ref_->target());
    shared_ptr<Dvec> rdm2d = r->rdm2deriv(ref_->target());
    shared_ptr<Dvec> rdm3d, rdm4d;
    // RDM4 is contracted a priori by the Fock operator
    tie(rdm3d, rdm4d) = r->rdm34deriv(ref_->target(), fockact);
    assert(rdm3d->ij() == rdm4d->ij());

    vector<IndexRange> o3 = {ci_, active_, active_};
    vector<IndexRange> o5 = {ci_, active_, active_, active_, active_};
    vector<IndexRange> o7 = {ci_, active_, active_, active_, active_, active_, active_};
    rdm1deriv_ = make_shared<Tensor>(o3);
    rdm2deriv_ = make_shared<Tensor>(o5);
    rdm3deriv_ = make_shared<Tensor>(o7);
    rdm4deriv_ = make_shared<Tensor>(o7);

    const int nclo = ref_->nclosed();
    vector<int> inpoff3(2,nclo); inpoff3.push_back(0);
    vector<int> inpoff5(4,nclo); inpoff5.push_back(0);
    vector<int> inpoff7(6,nclo); inpoff7.push_back(0);

    const btas::CRange<3> range3(rdm1d->extent(0)*rdm1d->extent(1), nact, nact);
    const btas::CRange<5> range5(rdm2d->extent(0)*rdm2d->extent(1), nact, nact, nact, nact);
    const btas::CRange<7> range7(rdm3d->extent(0)*rdm3d->extent(1), nact, nact, nact, nact, nact, nact);

    rdm1d->resize(range3);
    rdm2d->resize(range5);
    rdm3d->resize(range7);
    rdm4d->resize(range7);
    fill_block<3>(rdm1deriv_, rdm1d, inpoff3, vector<IndexRange>(o3.rbegin(), o3.rend())); 
    fill_block<5>(rdm2deriv_, rdm2d, inpoff5, vector<IndexRange>(o5.rbegin(), o5.rend())); 
    fill_block<7>(rdm3deriv_, rdm3d, inpoff7, vector<IndexRange>(o7.rbegin(), o7.rend())); 
    fill_block<7>(rdm4deriv_, rdm4d, inpoff7, vector<IndexRange>(o7.rbegin(), o7.rend())); 

    timer.tick_print("RDM derivative evaluation");
  }


  // rdms.
  if (ref_->ciwfn()) {
    const int nclo = ref_->nclosed();

    rdm1_ = make_shared<Tensor>(vector<IndexRange>(2,active_));
    rdm2_ = make_shared<Tensor>(vector<IndexRange>(4,active_));
    rdm3_ = make_shared<Tensor>(vector<IndexRange>(6,active_));
    rdm4_ = make_shared<Tensor>(vector<IndexRange>(8,active_));

    shared_ptr<const RDM<1>> rdm1 = ref_->rdm1(ref_->target());
    shared_ptr<const RDM<2>> rdm2 = ref_->rdm2(ref_->target());
    shared_ptr<RDM<3>> rdm3;
    shared_ptr<RDM<4>> rdm4;
    tie(rdm3, rdm4) = ref_->compute_rdm34(ref_->target());

    fill_block<2>(rdm1_, rdm1, vector<int>(2,nclo), vector<IndexRange>(2,active_)); 
    fill_block<4>(rdm2_, rdm2, vector<int>(4,nclo), vector<IndexRange>(4,active_)); 
    fill_block<6>(rdm3_, rdm3, vector<int>(6,nclo), vector<IndexRange>(6,active_)); 
    fill_block<8>(rdm4_, rdm4, vector<int>(8,nclo), vector<IndexRange>(8,active_)); 

    timer.tick_print("RDM evaluation");

    // TODO this part has to be in the state loop
    // construct denominator
    denom_.push_back(make_shared<Denom>(*rdm1, *rdm2, *rdm3, *rdm4, *fockact));

    timer.tick_print("Denominator evaluation");
  }

  // set e0
  e0_ = compute_e0();
}


void SpinFreeMethod::print_iteration() const {
  cout << "      ---- iteration ----" << endl << endl;
  time_ = chrono::high_resolution_clock::now();
}


void SpinFreeMethod::print_iteration(const int i, const double en, const double err, const int ist) const {
  auto end = chrono::high_resolution_clock::now();
  const double tim = chrono::duration_cast<chrono::milliseconds>(end-time_).count() * 0.001;
  cout << "     " << setw(4) << i;
  if (ist >= 0)
    cout << setw(4) << ist;
  cout << setw(15) << fixed << setprecision(8) << en << setw(15) << fixed << setprecision(8) << err
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


shared_ptr<Tensor> SpinFreeMethod::init_amplitude() const {
  shared_ptr<Tensor> out = v2_->clone();
  auto put = [this, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = v2_->get_size_alloc(i0, i1, i2, i3);
    unique_ptr<double[]> buf(new double[size]);
    fill_n(buf.get(), size, 0.0);
    out->put_block(buf, i0, i1, i2, i3);
  };
  for (auto& i3 : virt_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          put(i0, i1, i2, i3);
  for (auto& i2 : active_)
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          put(i0, i1, i2, i3);
  for (auto& i0 : active_)
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          put(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          put(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          put(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_) {
          put(i0, i1, i2, i3);
          put(i2, i1, i0, i3);
        }
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i1 : virt_)
          put(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i0 : active_)
        for (auto& i2 : closed_)
          put(i0, i1, i2, i3);
  return out;
}


double SpinFreeMethod::dot_product_transpose(shared_ptr<const MultiTensor> r, shared_ptr<const MultiTensor> t2) const {
  assert(r->nref() == t2->nref());
  double out = 0.0;
  for (int i = 0; i != r->nref(); ++i)
    out += r->fac(i) * t2->fac(i);
  for (auto& i : *r)
    for (auto& j : *t2)
      out += dot_product_transpose(i, j);
  return out;
}


double SpinFreeMethod::dot_product_transpose(shared_ptr<const Tensor> r, shared_ptr<const Tensor> t2) const {
  auto prod = [this, &r, &t2](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = r->get_size_alloc(i2, i3, i0, i1);
    if (size == 0) return 0.0;

    unique_ptr<double[]> tmp0 = t2->get_block(i0, i1, i2, i3);
    unique_ptr<double[]> tmp1(new double[size]);
    sort_indices<2,3,0,1,0,1,1,1>(tmp0.get(), tmp1.get(), i0.size(), i1.size(), i2.size(), i3.size());

    return blas::dot_product(tmp1.get(), size, r->get_block(i2, i3, i0, i1).get());
  };
  double out = 0.0;
  for (auto& i3 : virt_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          out += prod(i0, i1, i2, i3);
  for (auto& i2 : active_)
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          out += prod(i0, i1, i2, i3);
  for (auto& i0 : active_)
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          out += prod(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          out += prod(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          out += prod(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_) {
          out += prod(i0, i1, i2, i3);
          out += prod(i2, i1, i0, i3);
        }
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i1 : virt_)
          out += prod(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i0 : active_)
        for (auto& i2 : closed_)
          out += prod(i0, i1, i2, i3);
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
