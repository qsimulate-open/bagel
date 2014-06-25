//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/util/taskqueue.h>
#include <src/df/dfblock_london.h>
#include <src/integral/libint/libint.h>
#include <src/integral/rys/eribatch.h>
#include <src/util/constants.h>
#include <src/util/simple.h>

using namespace bagel;
using namespace std;


// TODO this is a temp fix
shared_ptr<DFBlock_London> DFBlock_London::transform_second(std::shared_ptr<const ZMatView> cmat2, const bool trans) const {
  auto c2 = make_shared<ZMatrix>(*cmat2);
  return transform_second(c2, trans);
}


shared_ptr<DFBlock_London> DFBlock_London::transform_third(std::shared_ptr<const ZMatView> cmat2, const bool trans) const {
  auto c2 = make_shared<ZMatrix>(*cmat2);
  return transform_third(c2, trans);
}


shared_ptr<DFBlock_London> DFBlock_London::transform_second(std::shared_ptr<const ZMatrix> cmat2, const bool trans) const {
  /*****/
  // We need the conjugate of the coefficient here, because we're transforming the index in the bra.  This should not be needed for transform_third
  shared_ptr<const ZMatrix> intermediate = cmat2->transpose();
  shared_ptr<const ZMatrix> cmat = intermediate->transpose_conjg();
  /*****/

  assert(trans ? cmat->mdim() : cmat->ndim() == b1size());

  const complex<double>* const c = cmat->data();
  const int nocc = trans ? cmat->ndim() : cmat->mdim();

  // so far I only consider the following case
  assert(b1start_ == 0);
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize(), nocc, b2size(), astart_, 0, b2start_, averaged_);

  for (size_t i = 0; i != b2size(); ++i) {
    if (!trans)
      zgemm3m_("N", "N", asize(), nocc, b1size(), 1.0, data()+i*asize()*b1size(), asize(), c, b1size(), 0.0, out->data()+i*asize()*nocc, asize());
    else
      zgemm3m_("N", "C", asize(), nocc, b1size(), 1.0, data()+i*asize()*b1size(), asize(), c, nocc, 0.0, out->data()+i*asize()*nocc, asize());
  }
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::transform_third(std::shared_ptr<const ZMatrix> cmat, const bool trans) const {
  assert(trans ? cmat->mdim() : cmat->ndim() == b2size());
  const complex<double>* const c = cmat->data();
  const int nocc = trans ? cmat->ndim() : cmat->mdim();

  // so far I only consider the following case
  assert(b2start_ == 0);
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize(), b1size(), nocc, astart_, b1start_, 0, averaged_);

  if (!trans)
    zgemm3m_("N", "N", asize()*b1size(), nocc, b2size(), 1.0, data(), asize()*b1size(), c, b2size(), 0.0, out->data(), asize()*b1size());
  else  // trans -> back transform
    zgemm3m_("N", "C", asize()*b1size(), nocc, b2size(), 1.0, data(), asize()*b1size(), c, nocc, 0.0, out->data(), asize()*b1size());
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::clone() const {
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize(), b1size(), b2size(), astart_, b1start_, b2start_, averaged_);
  out->zero();
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::copy() const {
  return make_shared<DFBlock_London>(*this);
}


void DFBlock_London::add_direct_product(const shared_ptr<const ZMatrix> a, const shared_ptr<const ZMatrix> b, const double fac) {
  assert(asize() == a->ndim() && b1size()*b2size() == b->size());
  zgeru_(asize(), b1size()*b2size(), fac, a->data(), 1, b->data(), 1, data(), asize());
}


shared_ptr<DFBlock_London> DFBlock_London::swap() const {
  auto out = make_shared<DFBlock_London>(adist_shell_, adist_, asize(), b2size(), b1size(), astart_, b2start_, b1start_, averaged_);
  for (size_t b2 = b2start_; b2 != b2start_+b2size(); ++b2)
    for (size_t b1 = b1start_; b1 != b1start_+b1size(); ++b1)
      copy_n(data()+asize()*(b1+b1size()*b2), asize(), out->data()+asize()*(b2+b2size()*b1));
  return out;
}


/*
shared_ptr<DFBlock_London> DFBlock_London::apply_rhf_2RDM(const double scale_exch) const {
  assert(b1size() == b2size());
  const int nocc = b1size();
  shared_ptr<DFBlock_London> out = clone();
  out->zero();
  // exchange contributions
  out->ax_plus_y(-2.0*scale_exch, *this);
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(1.0, data()+asize()*(i+nocc*i), asize(), diagsum.get());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(4.0, diagsum.get(), asize(), out->get()+asize()*(i+nocc*i));
  return out;
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DFBlock_London> DFBlock_London::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(b1size() == b2size());
  const int nocc = b1size();
  shared_ptr<DFBlock_London> out = clone();
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    zgemm3m_("N", "N", asize()*nocc, nocc, nocc, 1.0, data(), asize()*nocc, amat, nocc, 0.0, d2.get(), asize()*nocc);
    for (int i = 0; i != nocc; ++i)
      zgemm3m_("N", "N", asize(), nocc, nocc, -1.0, d2.get()+asize()*nocc*i, asize(), amat, nocc, 0.0, out->get()+asize()*nocc*i, asize());
    zgemm3m_("N", "N", asize()*nocc, nocc, nocc, 1.0, data(), asize()*nocc, bmat, nocc, 0.0, d2.get(), asize()*nocc);
    for (int i = 0; i != nocc; ++i)
      zgemm3m_("N", "N", asize(), nocc, nocc, -1.0, d2.get()+asize()*nocc*i, asize(), bmat, nocc, 1.0, out->get()+asize()*nocc*i, asize());
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc];
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], data()+asize()*(i+nocc*i), asize(), diagsum.get());
  for (int i = 0; i != nocc; ++i)
    blas::ax_plus_y_n(sum[i], diagsum.get(), asize(), out->get()+asize()*(i+nocc*i));
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size() && b1size() == b2size());
  // checking if natural orbitals...
  bool natural = true;
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__*100)
      natural = false;
  }
  shared_ptr<DFBlock_London> out = clone();
  out->zero();
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      daxpy_(asize(), -2.0, data()+asize()*(j+b1size()*i), 1, out->get()+asize()*(j+b1size()*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[asize()]);
  fill_n(diagsum.get(), asize(), 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 1.0, data()+asize()*(i+b1size()*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 4.0, diagsum.get(), 1, out->get()+asize()*(i+b1size()*i), 1);

  // act-act part
  // compress
  unique_ptr<double[]> buf(new double[nact*nact*asize()]);
  unique_ptr<double[]> buf2(new double[nact*nact*asize()]);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(data()+asize()*(j+nclosed+b1size()*(i+nclosed)), asize(), buf.get()+asize()*(j+nact*i));
  // multiply
  zgemm3m_("N", "N", asize(), nact*nact, nact*nact, 1.0, buf.get(), asize(), rdm, nact*nact, 0.0, buf2.get(), asize());
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      copy_n(buf2.get()+asize()*(j+nact*i), asize(), out->get()+asize()*(j+nclosed+b1size()*(i+nclosed)));

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  if (natural) {
    for (int i = 0; i != nact; ++i)
      daxpy_(asize(), 2.0*rdm1[i+nact*i], diagsum.get(), 1, out->get()+asize()*(i+nclosed+b1size()*(i+nclosed)), 1);
  } else {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        daxpy_(asize(), 2.0*rdm1[j+nact*i], diagsum.get(), 1, out->get()+asize()*(j+nclosed+b1size()*(i+nclosed)), 1);
  }
  unique_ptr<double[]> diagsum2(new double[asize()]);
  zgemv_("N", asize(), nact*nact, 1.0, buf.get(), asize(), rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize(), 2.0, diagsum2.get(), 1, out->get()+asize()*(i+b1size()*i), 1);
  // exchange contribution
  if (natural) {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j) {
        daxpy_(asize(), -rdm1[i+nact*i], data()+asize()*(j+b1size()*(i+nclosed)), 1, out->get()+asize()*(j+b1size()*(i+nclosed)), 1);
        daxpy_(asize(), -rdm1[i+nact*i], data()+asize()*(i+nclosed+b1size()*j), 1, out->get()+asize()*(i+nclosed+b1size()*j), 1);
      }
  } else {
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i1 = 0; i1 != nact; ++i1)
        for (int j = 0; j != nclosed; ++j) { // TODO be careful when rdm1 is not symmetric (e.g., transition density matrices)
          daxpy_(asize(), -rdm1[i1+nact*i0], data()+asize()*(j+b1size()*(i0+nclosed)), 1, out->get()+asize()*(j+b1size()*(i1+nclosed)), 1);
          daxpy_(asize(), -rdm1[i1+nact*i0], data()+asize()*(i0+nclosed+b1size()*j), 1, out->get()+asize()*(i1+nclosed+b1size()*j), 1);
        }
  }
  return out;
}


shared_ptr<DFBlock_London> DFBlock_London::apply_2RDM(const double* rdm) const {
  shared_ptr<DFBlock_London> out = clone();
  zgemm3m_("N", "C", asize(), b1size()*b2size(), b1size()*b2size(), 1.0, data(), asize(), rdm, b1size()*b2size(), 0.0, out->get(), asize());
  return out;
}
*/


shared_ptr<ZMatrix> DFBlock_London::form_2index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (asize() != o->asize() || (b1size() != o->b1size() && b2size() != o->b2size())) throw logic_error("illegal call of DFBlock_London::form_2index");
  shared_ptr<ZMatrix> target;

  if (b1size() == o->b1size()) {
    target = make_shared<ZMatrix>(b2size(),o->b2size());
    zgemm3m_("C", "N", b2size(), o->b2size(), asize()*b1size(), a, data(), asize()*b1size(), o->data(), asize()*b1size(), 0.0, target->data(), b2size());
  } else {
    assert(b2size() == o->b2size());
    target = make_shared<ZMatrix>(b1size(),o->b1size());
    for (int i = 0; i != b2size(); ++i)
      zgemm3m_("C", "N", b1size(), o->b1size(), asize(), a, data()+i*asize()*b1size(), asize(), o->data()+i*asize()*o->b1size(), asize(), 1.0, target->data(), b1size());
  }

  return target;
}


shared_ptr<ZMatrix> DFBlock_London::form_4index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (asize() != o->asize()) throw logic_error("illegal call of DFBlock_London::form_4index");
  auto target = make_shared<ZMatrix>(b1size()*b2size(), o->b1size()*o->b2size());
  zgemm3m_("C", "N", b1size()*b2size(), o->b1size()*o->b2size(), asize(), a, data(), asize(), o->data(), asize(), 0.0, target->data(), b1size()*b2size());
  return target;
}


// slowest index of o is fixed to n
shared_ptr<ZMatrix> DFBlock_London::form_4index_1fixed(const shared_ptr<const DFBlock_London> o, const double a, const size_t n) const {
  if (asize() != o->asize()) throw logic_error("illegal call of DFBlock_London::form_4index_1fixed");
  auto target = make_shared<ZMatrix>(b2size()*b1size(), o->b1size());
  zgemm3m_("C", "N", b1size()*b2size(), o->b1size(), asize(), a, data(), asize(), o->data()+n*asize()*o->b1size(), asize(), 0.0, target->data(), b1size()*b2size());
  return target;
}


shared_ptr<ZMatrix> DFBlock_London::form_aux_2index(const shared_ptr<const DFBlock_London> o, const double a) const {
  if (b1size() != o->b1size() || b2size() != o->b2size()) throw logic_error("illegal call of DFBlock_London::form_aux_2index");
  auto target = make_shared<ZMatrix>(asize(), o->asize());
  zgemm3m_("N", "C", asize(), o->asize(), b1size()*b2size(), a, data(), asize(), o->data(), o->asize(), 0.0, target->data(), asize());
  return target;
}


shared_ptr<ZVectorB> DFBlock_London::form_vec(const shared_ptr<const ZMatrix> den) const {
  auto out = make_shared<ZVectorB>(asize());
  auto dfv = btas::group(*this,  1, 3);
  auto denv = btas::group(*den, 0, 2);
  btas::contract(1.0, dfv, {0,1}, denv, {1}, 0.0, *out, {0});
  return out;
}


shared_ptr<ZMatrix> DFBlock_London::form_mat(const complex<double>* fit) const {
  auto out = make_shared<ZMatrix>(b1size(),b2size());
  zgemv_("T", asize(), b1size()*b2size(), 1.0, data(), asize(), fit, 1, 0.0, out->data(), 1);
  return out;
}


void DFBlock_London::contrib_apply_J(const shared_ptr<const DFBlock_London> o, const shared_ptr<const ZMatrix> d) {
  if (b1size() != o->b1size() || b2size() != o->b2size()) throw logic_error("illegal call of DFBlock_London::contrib_apply_J");
  zgemm3m_("N", "N", asize(), b1size()*b2size(), o->asize(), 1.0, d->element_ptr(astart_, o->astart_), d->ndim(), o->data(), o->asize(),
                                                         1.0, data(), asize());
}


shared_ptr<ZMatrix> DFBlock_London::form_Dj(const shared_ptr<const ZMatrix> o, const int jdim) const {
  assert(o->size() == b1size()*b2size()*jdim);
  auto out = make_shared<ZMatrix>(asize(), jdim);
  zgemm3m_("N", "N", asize(), jdim, b1size()*b2size(), 1.0, data(), asize(), o->data(), b1size()*b2size(), 0.0, out->data(), asize());
  return out;
}


shared_ptr<btas::Tensor3<std::complex<double>>> DFBlock_London::get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize() || jfen > b1size() || kfen > b2size())
    throw logic_error("illegal call of DFBlock_London::get_block");

  // TODO we need 3-index tensor class here!
  auto out = make_shared<btas::Tensor3<std::complex<double>>>(i, j*k);
  for (int kk = ksta; kk != kfen; ++kk)
    for (int jj = jsta; jj != jfen; ++jj)
      for (int ii = ista; ii != ifen; ++ii)
        (*out)(ii-ista, jj-jsta, kk-ksta) = (*this)(ii, jj, kk);

  return out;
}


shared_ptr<btas::Tensor3<std::complex<double>>> DFBlock_London::get_block_conj(const int ist, const int i, const int jst, const int j, const int kst, const int k) const {
  if (ist != 0 || jst != 0 || kst != 0 || astart_ != 0 || b1start_ != 0 || b2start_ != 0) throw logic_error("DFBlock_London::get_block_conj currently is not designed to work with >1 block");
  if (b1size() != b2size()) throw logic_error ("DFBLock_London::get_block_conj assumes b1 and b2 contain the same set of basis functions");

  const int ista = ist - astart_;
  const int jsta = jst - b1start_;
  const int ksta = kst - b2start_;
  const int ifen = ist + i - astart_;
  const int jfen = jst + j - b1start_;
  const int kfen = kst + k - b2start_;
  if (ista < 0 || jsta < 0 || ksta < 0 || ifen > asize() || jfen > b1size() || kfen > b2size())
    throw logic_error("illegal call of DFBlock_London::get_block");

  // TODO we need 3-index tensor class here!
  auto out = make_shared<btas::Tensor3<std::complex<double>>>(i, j, k);

  for (int kk = ksta; kk != kfen; kk++)
    for (int jj = jsta; jj != jfen; jj++)
      for (int ii = ista; ii != ifen; ii++)
        (*out)(ii-ista, jj-jsta, kk-ksta) = conj((*this)(ii, kk, jj));

  return out;
}
