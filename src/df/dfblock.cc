//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/df/dfblock.h>

using namespace bagel;
using namespace std;
using namespace btas;

shared_ptr<DFBlock> DFBlock::transform_second(std::shared_ptr<const TensorView2<double>> cmat, const bool trans) const {
  assert(trans ? cmat->extent(1) : cmat->extent(0) == b1size());
  assert(cmat->range().ordinal().contiguous());

  // so far I only consider the following case
  assert(b1start_ == 0);
  const int nocc = trans ? cmat->extent(0) : cmat->extent(1);
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), nocc, b2size(), astart_, 0, b2start_, averaged_);

  if (!trans)
    contract(1.0, *this, {0,3,2}, *cmat, {3,1}, 0.0, *out, {0,1,2});
  else
    contract(1.0, *this, {0,3,2}, *cmat, {1,3}, 0.0, *out, {0,1,2});

  return out;
}


shared_ptr<DFBlock> DFBlock::transform_third(std::shared_ptr<const TensorView2<double>> cmat, const bool trans) const {
  assert(trans ? cmat->extent(1) : cmat->extent(0) == b2size());
  assert(cmat->range().ordinal().contiguous());

  // so far I only consider the following case
  assert(b2start_ == 0);
  const int nocc = trans ? cmat->extent(0) : cmat->extent(1);
  auto out = make_shared<DFBlock>(adist_shell_, adist_, asize(), b1size(), nocc, astart_, b1start_, 0, averaged_);

  if (!trans)
    contract(1.0, *this, {0,1,3}, *cmat, {3,2}, 0.0, *out, {0,1,2});
  else  // trans -> back transform
    contract(1.0, *this, {0,1,3}, *cmat, {2,3}, 0.0, *out, {0,1,2});

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


void DFBlock::add_direct_product(const shared_ptr<const VectorB> a, const shared_ptr<const Matrix> b, const double fac) {
  assert(asize() == a->size() && b1size()*b2size() == b->size());
  dger_(asize(), b1size()*b2size(), fac, a->data(), 1, b->data(), 1, data(), asize());
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
      daxpy_(asize(), -2.0, data()+asize()*(j+b1size()*i), 1, out->data()+asize()*(j+b1size()*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 1.0, data()+asize()*(i+b1size()*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 4.0, diagsum.get(), 1, out->data()+asize()*(i+b1size()*i), 1);

  // act-act part
  // compress
  auto low = {0, nclosed, nclosed};
  auto up = {static_cast<int>(asize()), nclosed+nact, nclosed+nact};
  Tensor3<double> buf = TensorView3<double>(range().slice(low, up), storage());
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
      daxpy_(asize(), 2.0*rdm1(i, i), diagsum.get(), 1, out->data()+asize()*(i+nclosed+b1size()*(i+nclosed)), 1);
  } else {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        daxpy_(asize(), 2.0*rdm1(j, i), diagsum.get(), 1, out->data()+asize()*(j+nclosed+b1size()*(i+nclosed)), 1);
  }
  VectorB diagsum2(asize());
  contract(1.0, group(buf,1,3), {0,1}, group(rdm1,0,2), {1}, 0.0, diagsum2, {0});

  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 2.0, diagsum2.data(), 1, out->data()+asize()*(i+b1size()*i), 1);
  // exchange contribution
  if (natural) {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
        daxpy_(asize(), -rdm1(i, i), data()+asize()*(j+b1size()*(i+nclosed)), 1, out->data()+asize()*(j+b1size()*(i+nclosed)), 1);
        daxpy_(asize(), -rdm1(i, i), data()+asize()*(i+nclosed+b1size()*j), 1, out->data()+asize()*(i+nclosed+b1size()*j), 1);
      }
  } else {
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i1 = 0; i1 != nact; ++i1)
        for (int j = 0; j != nclosed; ++j) { // TODO be careful when rdm1 is not symmetric (e.g., transition density matrices)
          daxpy_(asize(), -rdm1(i1, i0), data()+asize()*(j+b1size()*(i0+nclosed)), 1, out->data()+asize()*(j+b1size()*(i1+nclosed)), 1);
          daxpy_(asize(), -rdm1(i1, i0), data()+asize()*(i0+nclosed+b1size()*j), 1, out->data()+asize()*(i1+nclosed+b1size()*j), 1);
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

  // TODO we need 3-index tensor class here!
  auto out = make_shared<Tensor3<double>>(i, j, k);
  for (int kk = ksta; kk != kfen; ++kk)
    for (int jj = jsta; jj != jfen; ++jj)
      for (int ii = ista; ii != ifen; ++ii)
        (*out)(ii-ista, jj-jsta, kk-ksta) = (*this)(ii, jj, kk);

  return out;
}
