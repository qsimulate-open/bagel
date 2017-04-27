//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfblock.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/df/dfblock.h>

using namespace bagel;
using namespace std;
using namespace btas;

// construction of a block from AO integrals
DFBlock::DFBlock(shared_ptr<const StaticDist> adist_shell, shared_ptr<const StaticDist> adist,
             const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s, const bool averaged)
 : btas::Tensor3<double>(max(adist_shell->size(mpi__->rank()), max(adist->size(mpi__->rank()), a)), b1, b2),
   adist_shell_(adist_shell), adist_(adist), averaged_(averaged), astart_(as), b1start_(b1s), b2start_(b2s) {

  assert(asize() == adist_shell->size(mpi__->rank()) || asize() == adist_->size(mpi__->rank()) || asize() == adist_->nele());

  // resize to the current size (moving the end pointer)
  const btas::CRange<3> range(a, b1, b2);
  this->resize(range);
}


DFBlock::DFBlock(const DFBlock& o)
 : btas::Tensor3<double>(max(o.adist_shell_->size(mpi__->rank()), max(o.adist_->size(mpi__->rank()), o.asize())), o.b1size(), o.b2size()),
   adist_shell_(o.adist_shell_), adist_(o.adist_), averaged_(o.averaged_), astart_(o.astart_), b1start_(o.b1start_), b2start_(o.b2start_) {

  // resize to the current size
  const btas::CRange<3> range(o.asize(), o.b1size(), o.b2size());
  this->resize(range);

  btas::Tensor3<double>::operator=(o);
}


shared_ptr<DFBlock> DFBlock::transform_second(const MatView cmat, const bool trans) const {
  assert(trans ? cmat.extent(1) : cmat.extent(0) == b1size());
  assert(cmat.range().ordinal().contiguous());

  // so far I only consider the following case
  assert(b1start_ == 0);
  const int nocc = trans ? cmat.extent(0) : cmat.extent(1);
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), nocc, b2size(), astart_, 0, b2start_, averaged_);

  if (!trans)
    contract(1.0, *this, {0,3,2}, cmat, {3,1}, 0.0, *out, {0,1,2});
  else
    contract(1.0, *this, {0,3,2}, cmat, {1,3}, 0.0, *out, {0,1,2});

  return out;
}


shared_ptr<DFBlock> DFBlock::transform_third(const MatView cmat, const bool trans) const {
  assert(trans ? cmat.extent(1) : cmat.extent(0) == b2size());
  assert(cmat.range().ordinal().contiguous());

  // so far I only consider the following case
  assert(b2start_ == 0);
  const int nocc = trans ? cmat.extent(0) : cmat.extent(1);
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), b1size(), nocc, astart_, b1start_, 0, averaged_);

  if (!trans)
    contract(1.0, *this, {0,1,3}, cmat, {3,2}, 0.0, *out, {0,1,2});
  else  // trans -> back transform
    contract(1.0, *this, {0,1,3}, cmat, {2,3}, 0.0, *out, {0,1,2});

  return out;
}


shared_ptr<DFBlock> DFBlock::merge_b1(shared_ptr<const DFBlock> o) const {
  assert(asize() == o->asize() && b2size() == o->b2size());
  assert(astart() == o->astart() && b1start() == o->b1start() && b2start() == o->b2start());
  assert(adist_shell_ == o->adist_shell_ && adist_ == o->adist_ && averaged_ == o->averaged_);

  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), b1size() + o->b1size(), b2size(), astart_, b1start_, b2start_, averaged_);
  const int size1 = asize() * b1size();
  const int size2 = o->asize() * o->b1size();
  const int size3 = size1 + size2;
  for (int i = 0; i != b2size(); ++i) {
    copy_n(   data() + i*size1, size1, out->data() + i*size3);
    copy_n(o->data() + i*size2, size2, out->data() + i*size3 + size1);
  }
  return out;
}


shared_ptr<DFBlock> DFBlock::slice_b1(const int slice_start, const int slice_size) const {
  assert(slice_start >= 0 && slice_start + slice_size <= b1size());
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), slice_size, b2size(), astart_, b1start_, b2start_, averaged_);
  const int size1 = asize() * slice_size;
  const int size2 = asize() * b1size();
  for (int i = 0; i != b2size(); ++i) {
    copy_n(data() + asize()*slice_start + i*size2, size1, out->data() + i*size1);
  }
  return out;
}


shared_ptr<DFBlock> DFBlock::clone() const {
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), b1size(), b2size(), astart_, b1start_, b2start_, averaged_);
  out->zero();
  return out;
}


shared_ptr<DFBlock> DFBlock::copy() const {
  return make_shared<DFBlock>(*this);
}


void DFBlock::add_direct_product(const VecView a, const MatView b, const double fac) {
  assert(asize() == a.size() && b1size()*b2size() == b.size());
  dger_(asize(), b1size()*b2size(), fac, a.data(), 1, b.data(), 1, data(), asize());
}


shared_ptr<DFBlock> DFBlock::swap() const {
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), b2size(), b1size(), astart_, b2start_, b1start_, averaged_);
  for (size_t b2 = b2start_; b2 != b2start_+b2size(); ++b2)
    for (size_t b1 = b1start_; b1 != b1start_+b1size(); ++b1)
      copy_n(data()+asize()*(b1+b1size()*b2), asize(), out->data()+asize()*(b2+b2size()*b1));
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_rhf_2RDM(const double scale_exch) const {
  assert(b1size() == b2size());
  const int nocc = b1size();
  shared_ptr<DFBlock> out = clone();
  out->zero();
  // exchange contributions
  out->ax_plus_y(-2.0*scale_exch, *this);
  // coulomb contributions (diagonal to diagonal)
  VectorB diagsum(asize());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(1.0, data()+asize()*(i+nocc*i), asize(), diagsum.data());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(4.0, diagsum.data(), asize(), out->data()+asize()*(i+nocc*i));
  return out;
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DFBlock> DFBlock::apply_uhf_2RDM(const Tensor2<double>& amat, const Tensor2<double>& bmat) const {
  assert(b1size() == b2size());
  const int nocc = b1size();
  shared_ptr<DFBlock> out = clone();
  {
    auto d2 = clone();
    // exchange contributions
    contract( 1.0, *this, {0,1,2}, amat, {2,3}, 0.0,  *d2, {0,1,3});
    contract(-1.0,   *d2, {0,1,2}, amat, {1,3}, 0.0, *out, {0,3,2});
    contract( 1.0, *this, {0,1,2}, bmat, {2,3}, 0.0,  *d2, {0,1,3});
    contract(-1.0,   *d2, {0,1,2}, bmat, {1,3}, 1.0, *out, {0,3,2});
  }

  VectorB sum(nocc);
  for (int i = 0; i != nocc; ++i)
    sum[i] = amat(i,i) + bmat(i,i);

  // coulomb contributions (diagonal to diagonal)
  VectorB diagsum(asize());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], data()+asize()*(i+nocc*i), asize(), diagsum.data());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], diagsum.data(), asize(), out->data()+asize()*(i+nocc*i));
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_2RDM(const Tensor4<double>& rdm, const Tensor2<double>& rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size() && b1size() == b2size());
  // checking if natural orbitals...
  bool natural = true;
  {
    const double a = ddot_(nact*nact, rdm1.data(), 1, rdm1.data(), 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1(i,i)*rdm1(i,i);
    if (fabs(a-sum) > numerical_zero__*100)
      natural = false;
  }
  shared_ptr<DFBlock> out = clone();
  out->zero();
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      blas::ax_plus_y_n(-2.0, data()+asize()*(j+b1size()*i), asize(), out->data()+asize()*(j+b1size()*i));
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nclosed; ++i)
    blas::ax_plus_y_n(1.0, data()+asize()*(i+b1size()*i), asize(), diagsum.get());
  for (int i = 0; i != nclosed; ++i)
    blas::ax_plus_y_n(4.0, diagsum.get(), asize(), out->data()+asize()*(i+b1size()*i));

  // act-act part
  // compress
  auto low = {0, nclosed, nclosed};
  auto up = {static_cast<int>(asize()), nclosed+nact, nclosed+nact};
  Tensor3<double> buf = make_view(range().slice(low, up), storage());
  Tensor3<double> buf2(asize(), nact, nact);
  {
    auto rdm2v = group(group(rdm,2,4),0,2);
    auto buf2v = group(buf2,1,3);
    contract(1.0, group(buf,1,3), {0,1}, rdm2v, {1,2}, 0.0, buf2v, {0,2});
  }
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(buf2.data()+asize()*(j+nact*i), asize(), out->data()+asize()*(j+nclosed+b1size()*(i+nclosed)));

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  if (natural) {
    for (int i = 0; i != nact; ++i)
      blas::ax_plus_y_n(2.0*rdm1(i, i), diagsum.get(), asize(), out->data()+asize()*(i+nclosed+b1size()*(i+nclosed)));
  } else {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        blas::ax_plus_y_n(2.0*rdm1(j, i), diagsum.get(), asize(), out->data()+asize()*(j+nclosed+b1size()*(i+nclosed)));
  }
  VectorB diagsum2(asize());
  contract(1.0, group(buf,1,3), {0,1}, group(rdm1,0,2), {1}, 0.0, diagsum2, {0});

  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 2.0, diagsum2.data(), 1, out->data()+asize()*(i+b1size()*i), 1);
  // exchange contribution
  if (natural) {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
        blas::ax_plus_y_n(-rdm1(i, i), data()+asize()*(j+b1size()*(i+nclosed)), asize(), out->data()+asize()*(j+b1size()*(i+nclosed)));
        blas::ax_plus_y_n(-rdm1(i, i), data()+asize()*(i+nclosed+b1size()*j), asize(), out->data()+asize()*(i+nclosed+b1size()*j));
      }
  } else {
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i1 = 0; i1 != nact; ++i1)
        for (int j = 0; j != nclosed; ++j) { // TODO be careful when rdm1 is not symmetric (e.g., transition density matrices)
          blas::ax_plus_y_n(-rdm1(i1, i0), data()+asize()*(j+b1size()*(i0+nclosed)), asize(), out->data()+asize()*(j+b1size()*(i1+nclosed)));
          blas::ax_plus_y_n(-rdm1(i1, i0), data()+asize()*(i0+nclosed+b1size()*j), asize(), out->data()+asize()*(i1+nclosed+b1size()*j));
        }
  }
  return out;
}

shared_ptr<DFBlock> DFBlock::apply_2RDM_tr(const Tensor4<double>& rdm, const Tensor2<double>& rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size() && b1size() == b2size());
  // checking if natural orbitals...

  bool natural = true;
  {
    const double a = ddot_(nact*nact, rdm1.data(), 1, rdm1.data(), 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1(i,i)*rdm1(i,i);
    if (fabs(a-sum) > numerical_zero__*100)
      natural = false;
  }
  shared_ptr<DFBlock> out = clone();
  out->zero();

  // initialize diagonal elements
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nclosed; ++i)
    blas::ax_plus_y_n(1.0, data()+asize()*(i+b1size()*i), asize(), diagsum.get());

  // act-act part
  // compress
  auto low = {0, nclosed, nclosed};
  auto up = {static_cast<int>(asize()), nclosed+nact, nclosed+nact};
  Tensor3<double> buf = make_view(range().slice(low, up), storage());
  Tensor3<double> buf2(asize(), nact, nact);
  {
    auto rdm2v = group(group(rdm,2,4),0,2);
    auto buf2v = group(buf2,1,3);
    contract(1.0, group(buf,1,3), {0,1}, rdm2v, {1,2}, 0.0, buf2v, {0,2});
  }
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(buf2.data()+asize()*(j+nact*i), asize(), out->data()+asize()*(j+nclosed+b1size()*(i+nclosed)));

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  if (natural) {
    for (int i = 0; i != nact; ++i)
      blas::ax_plus_y_n(2.0*rdm1(i, i), diagsum.get(), asize(), out->data()+asize()*(i+nclosed+b1size()*(i+nclosed)));
  } else {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j) {
        blas::ax_plus_y_n(2.0*rdm1(j, i), diagsum.get(), asize(), out->data()+asize()*(j+nclosed+b1size()*(i+nclosed)));
      }
  }

  VectorB diagsum2(asize());
  contract(1.0, group(buf,1,3), {0,1}, group(rdm1,0,2), {1}, 0.0, diagsum2, {0});

  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 2.0, diagsum2.data(), 1, out->data()+asize()*(i+b1size()*i), 1);
  // exchange contribution
  if (natural) {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
        blas::ax_plus_y_n(-rdm1(i, i), data()+asize()*(j+b1size()*(i+nclosed)), asize(), out->data()+asize()*(j+b1size()*(i+nclosed)));
        blas::ax_plus_y_n(-rdm1(i, i), data()+asize()*(i+nclosed+b1size()*j), asize(), out->data()+asize()*(i+nclosed+b1size()*j));
      }
  } else {
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i1 = 0; i1 != nact; ++i1) {
        double gamma = -(rdm1(i1, i0) + rdm1(i0, i1)) * .5;
        for (int j = 0; j != nclosed; ++j) {
          blas::ax_plus_y_n(gamma, data()+asize()*(i1+nclosed+b1size()*j), asize(), out->data()+asize()*(j+b1size()*(i0+nclosed)));
          blas::ax_plus_y_n(gamma, data()+asize()*(j+b1size()*(i1+nclosed)), asize(), out->data()+asize()*(i0+nclosed+b1size()*j));
        }
      }
  }
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_2RDM(const Tensor4<double>& rdm) const {
  shared_ptr<DFBlock> out = clone();
  auto rdmv = group(group(rdm,  2, 4), 0, 2);
  auto outv = group(*out,1,3);
  contract(1.0, group(*this,1,3), {0,2}, rdmv, {2,1}, 0.0, outv, {0,1});
  return out;
}


shared_ptr<Matrix> DFBlock::form_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize() != o->asize() || (b1size() != o->b1size() && b2size() != o->b2size())) throw logic_error("illegal call of DFBlock::form_2index");
  shared_ptr<Matrix> target;

  if (b1size() == o->b1size()) {
    target = make_shared<Matrix>(b2size(),o->b2size());
    contract(a, *this, {2,3,0}, *o, {2,3,1}, 0.0, *target, {0,1});
  } else {
    assert(b2size() == o->b2size());
    target = make_shared<Matrix>(b1size(),o->b1size());
    contract(a, *this, {2,0,3}, *o, {2,1,3}, 0.0, *target, {0,1});
  }

  return target;
}


shared_ptr<Matrix> DFBlock::form_4index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize() != o->asize()) throw logic_error("illegal call of DFBlock::form_4index");
  auto target = make_shared<Matrix>(b1size()*b2size(), o->b1size()*o->b2size());
  contract(a, group(*this,1,3), {1,0}, group(*o,1,3), {1,2}, 0.0, *target, {0,2});
  return target;
}


// slowest index of o is fixed to n
shared_ptr<Matrix> DFBlock::form_4index_1fixed(const shared_ptr<const DFBlock> o, const double a, const size_t n) const {
  if (asize() != o->asize()) throw logic_error("illegal call of DFBlock::form_4index_1fixed");
  auto target = make_shared<Matrix>(b2size()*b1size(), o->b1size());
  dgemm_("T", "N", b1size()*b2size(), o->b1size(), asize(), a, data(), asize(), o->data()+n*asize()*o->b1size(), asize(), 0.0, target->data(), b1size()*b2size());
  return target;
}


shared_ptr<Matrix> DFBlock::form_4index_diagonal() const {
  auto target = make_shared<Matrix>(b1size(), b2size());
  for (int i = 0; i != b2size(); ++i)
    for (int j = 0; j != b1size(); ++j)
      target->element(j, i) = blas::dot_product(data()+asize()*(j+b1size()*i), asize(), data()+asize()*(j+b1size()*i));
  return target;
}


shared_ptr<Matrix> DFBlock::form_4index_diagonal_part() const {
  // not very good code (assuming that b1size() is small)
  auto target = make_shared<Matrix>(b1size()*b1size(), b2size());
  for (int i = 0; i != b2size(); ++i)
    for (int j = 0; j != b1size(); ++j)
      for (int k = 0; k != b1size(); ++k)
         target->element(k+b1size()*j, i) = blas::dot_product(data()+asize()*(k+b1size()*i), asize(), data()+asize()*(j+b1size()*i));
  return target;
}


shared_ptr<Matrix> DFBlock::form_aux_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (b1size() != o->b1size() || b2size() != o->b2size()) throw logic_error("illegal call of DFBlock::form_aux_2index");
  auto target = make_shared<Matrix>(asize(), o->asize());
  contract(a, *this, {0,2,3}, *o, {1,2,3}, 0.0, *target, {0,1});
  return target;
}


shared_ptr<VectorB> DFBlock::form_vec(const shared_ptr<const Matrix> den) const {
  auto out = make_shared<VectorB>(asize());
  contract(1.0, group(*this,1,3), {0,1}, group(*den,0,2), {1}, 0.0, *out, {0});
  return out;
}


shared_ptr<Matrix> DFBlock::form_mat(const Tensor1<double>& fit) const {
  auto out = make_shared<Matrix>(b1size(), b2size());
  auto outv = group(*out,0,2);
  contract(1.0, group(*this,1,3), {1,0}, fit, {1}, 0.0, outv, {0});
  return out;
}


void DFBlock::contrib_apply_J(const shared_ptr<const DFBlock> o, const shared_ptr<const Matrix> d) {
  if (b1size() != o->b1size() || b2size() != o->b2size()) throw logic_error("illegal call of DFBlock::contrib_apply_J");
  assert(astart_ == 0 && o->astart_ == 0);
  contract(1.0, *d, {0,3}, *o, {3,1,2}, 1.0, *this, {0,1,2});
}


shared_ptr<Matrix> DFBlock::form_Dj(const shared_ptr<const Matrix> o, const int jdim) const {
  assert(o->size() == b1size()*b2size()*jdim);
  auto out = make_shared<Matrix>(asize(), jdim);
  contract(1.0, group(*this,1,3), {0,1}, *o, {1,2}, 0.0, *out, {0,2});
  return out;
}


shared_ptr<Tensor3<double>> DFBlock::get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize() || jfen > b1size() || kfen > b2size())
    throw logic_error("illegal call of DFBlock::get_block");

  auto out = make_shared<Tensor3<double>>(i, j, k);
  for (int kk = ksta; kk != kfen; ++kk)
    for (int jj = jsta; jj != jfen; ++jj)
      copy_n(&((*this)(ista, jj, kk)), ifen-ista, &((*out)(0, jj-jsta, kk-ksta)));

  return out;
}


DFBlock& DFBlock::operator=(const DFBlock& o) {
  btas::Tensor3<double>::operator=(o);
  adist_shell_ = o.adist_shell_;
  adist_ = o.adist_;
  averaged_ = o.averaged_;
  astart_ = o.astart_;
  b1start_ = o.b1start_;
  b2start_ = o.b2start_;
  return *this;
}


DFBlock& DFBlock::operator=(DFBlock&& o) {
  btas::Tensor3<double>::operator=(move(o));
  adist_shell_ = o.adist_shell_;
  adist_ = o.adist_;
  averaged_ = o.averaged_;
  astart_ = o.astart_;
  b1start_ = o.b1start_;
  b2start_ = o.b2start_;
  return *this;
}


void DFBlock::symmetrize() {
  if (b1size() != b2size()) throw logic_error("illegal call of DFBlock::symmetrize()");
  const int n = b1size();
  for (int i = 0; i != n; ++i)
    for (int j = i; j != n; ++j) {
      blas::ax_plus_y_n(1.0, data()+asize()*(j+n*i), asize(), data()+asize()*(i+n*j));
      copy_n(data()+asize()*(i+n*j), asize(), data()+asize()*(j+n*i));
    }
}


void DFBlock::copy_block(shared_ptr<MatView> o, const int jdim, const size_t offset) {
  assert(o->size() == asize()*jdim);
  copy_n(o->data(), asize()*jdim, data()+offset);
}


void DFBlock::copy_block(MatView o, const int jdim, const size_t offset) {
  assert(o.size() == asize()*jdim);
  copy_n(o.data(), asize()*jdim, data()+offset);
}


void DFBlock::add_block(shared_ptr<MatView> o, const int jdim, const size_t offset, const double fac) {
  assert(o->size() == asize()*jdim);
  blas::ax_plus_y_n(fac, o->data(), asize()*jdim, data()+offset);
}


void DFBlock::add_block(MatView o, const int jdim, const size_t offset, const double fac) {
  assert(o.size() == asize()*jdim);
  blas::ax_plus_y_n(fac, o.data(), asize()*jdim, data()+offset);
}


// average the asize between MPI processes (block will be described by dist_)
void DFBlock::average() {
  if (averaged_) return;
  averaged_ = true;

  // first make a send and receive buffer
  const size_t o_start = astart_;
  const size_t o_end   = o_start + asize();
  const int myrank = mpi__->rank();
  size_t t_start, t_end;
  tie(t_start, t_end) = adist_->range(myrank);

  assert(o_end >= t_end);
  assert(o_start >= t_start);

  // TODO so far I am not considering the cases when data must be sent to the next neighbor; CAUTION
  const size_t asendsize = o_end - t_end;
  const size_t arecvsize = o_start - t_start;

  assert(asendsize < t_end-t_start && arecvsize < t_end-t_start);

  unique_ptr<double[]> sendbuf;
  unique_ptr<double[]> recvbuf;
  int sendtag = 0;
  int recvtag = 0;

  if (asendsize) {
    TaskQueue<CopyBlockTask<double>> task(b2size());

    sendbuf = unique_ptr<double[]>(new double[asendsize*b1size()*b2size()]);
    const size_t retsize = asize() - asendsize;
    for (size_t b2 = 0; b2 != b2size(); ++b2)
      task.emplace_back(data()+retsize+asize()*b1size()*b2, asize(), sendbuf.get()+asendsize*b1size()*b2, asendsize, asendsize, b1size());

    task.compute();

    // send to the next node
    sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size()*b2size(), myrank+1, myrank);
  }

  if (arecvsize) {
    recvbuf = unique_ptr<double[]>(new double[arecvsize*b1size()*b2size()]);
    // recv from the previous node
    recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size()*b2size(), myrank-1, myrank-1);
  }

  // second move local data
  if (arecvsize || asendsize) {
    const size_t t_size = t_end - t_start;
    const size_t retsize = asize() - asendsize;
    if (t_size <= asize()) {
      for (size_t i = 0; i != b1size()*b2size(); ++i) {
        if (i*asize() < (i+1)*t_size-retsize) {
          copy_backward(data()+i*asize(), data()+i*asize()+retsize, data()+(i+1)*t_size);
        } else if (i*asize() > (i+1)*t_size-retsize) {
          copy_n(data()+i*asize(), retsize, data()+(i+1)*t_size-retsize);
        }
      }
    } else {
      for (long long int i = b1size()*b2size()-1; i >= 0; --i) {
        assert(i*asize() < (i+1)*t_size-retsize);
        copy_backward(data()+i*asize(), data()+i*asize()+retsize, data()+(i+1)*t_size);
      }
    }
  }

  // set new astart_ and asize()
  astart_ = t_start;
  assert(this->storage().capacity() >= (t_end - t_start)*b1size()*b2size());
  const btas::CRange<3> range(t_end - t_start, b1size(), b2size());
  this->resize(range);

  // set received data
  if (arecvsize) {
    // wait for recv communication
    mpi__->wait(recvtag);

    TaskQueue<CopyBlockTask<double>> task(b2size());
    for (size_t b2 = 0; b2 != b2size(); ++b2)
      task.emplace_back(recvbuf.get()+arecvsize*b1size()*b2, arecvsize, data()+asize()*b1size()*b2, asize(), arecvsize, b1size());
    task.compute();
  }

  // wait for send communication
  if (asendsize) mpi__->wait(sendtag);
}


// reverse operation of average() function
void DFBlock::shell_boundary() {
  if (!averaged_) return;
  averaged_ = false;
  const size_t o_start = astart_;
  const size_t o_end = o_start + asize();
  const int myrank = mpi__->rank();
  size_t t_start, t_end;
  tie(t_start, t_end) = adist_shell_->range(myrank);

  const size_t asendsize = t_start - o_start;
  const size_t arecvsize = t_end - o_end;
  assert(t_start >= o_start && t_end >= o_end);

  unique_ptr<double[]> sendbuf, recvbuf;
  int sendtag = 0;
  int recvtag = 0;

  if (asendsize) {
    TaskQueue<CopyBlockTask<double>> task(b2size());
    sendbuf = unique_ptr<double[]>(new double[asendsize*b1size()*b2size()]);
    for (size_t b2 = 0; b2 != b2size(); ++b2)
      task.emplace_back(data()+asize()*b1size()*b2, asize(), sendbuf.get()+asendsize*b1size()*b2, asendsize, asendsize, b1size());

    task.compute();
    assert(myrank > 0);
    sendtag = mpi__->request_send(sendbuf.get(), asendsize*b1size()*b2size(), myrank-1, myrank);
  }
  if (arecvsize) {
    assert(myrank+1 < mpi__->size());
    recvbuf = unique_ptr<double[]>(new double[arecvsize*b1size()*b2size()]);
    recvtag = mpi__->request_recv(recvbuf.get(), arecvsize*b1size()*b2size(), myrank+1, myrank+1);
  }

  if (arecvsize || asendsize) {
    const size_t t_size = t_end - t_start;
    const size_t retsize = asize() - asendsize;
    assert(t_size >= retsize);
    if (t_size <= asize()) {
      for (size_t i = 0; i != b1size()*b2size(); ++i) {
        assert(i*asize()+asendsize > i*t_size);
        copy_n(data()+i*asize()+asendsize, retsize, data()+i*t_size);
      }
    } else {
      for (long long int i = b1size()*b2size()-1; i >= 0; --i) {
        if (i*asize()+asendsize > i*t_size) {
          copy_n(data()+i*asize()+asendsize, retsize, data()+i*t_size);
        } else if (i*asize()+asendsize < i*t_size) {
          copy_backward(data()+i*asize()+asendsize, data()+(i+1)*asize(), data()+i*t_size+retsize);
        }
      }
    }
  }

  // set new astart_ and asize()
  astart_ = t_start;
  assert(this->storage().capacity() >= (t_end - t_start)*b1size()*b2size());
  const btas::CRange<3> range(t_end - t_start, b1size(), b2size());
  this->resize(range);

  // set received data
  if (arecvsize) {
    // wait for recv communication
    mpi__->wait(recvtag);

    TaskQueue<CopyBlockTask<double>> task(b2size());
    for (size_t b2 = 0; b2 != b2size(); ++b2)
      task.emplace_back(recvbuf.get()+arecvsize*b1size()*b2, arecvsize, data()+asize()*b1size()*b2+(asize()-arecvsize), asize(), arecvsize, b1size());
    task.compute();
  }

  // wait for send communication
  if (asendsize) mpi__->wait(sendtag);
}
