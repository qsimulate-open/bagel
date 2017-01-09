//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldffull.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <src/df/reldffull.h>

using namespace std;
using namespace bagel;


RelDFFull::RelDFFull(shared_ptr<const RelDFHalf> df, shared_ptr<const ZMatrix> coeff) : RelDFBase(*df) {
  // Separate Coefficients into real and imaginary
  array<shared_ptr<const Matrix>,4> rcoeff;
  array<shared_ptr<const Matrix>,4> icoeff;
  assert(coeff->ndim() % 4 == 0);
  const size_t nbasis = coeff->ndim() / 4;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> oc = coeff->get_submatrix(i*nbasis, 0, nbasis, coeff->mdim());
    rcoeff[i] = oc->get_real_part();
    icoeff[i] = oc->get_imag_part();
  }
  init(df, rcoeff, icoeff);
}


RelDFFull::RelDFFull(shared_ptr<const RelDFHalf> df, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff) : RelDFBase(*df) {
  init(df, rcoeff, icoeff);
}


void RelDFFull::init(shared_ptr<const RelDFHalf> df, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff) {
  basis_ = df->basis();
  if (basis_.size() != 1)
    throw logic_error("RelDFFull should be called with basis_.size() == 1");

  const int index = basis_.front()->basis(1);

  // df->get_real() + df->get_imag() needed for 3-multiplication algorithm
  auto dfri = make_shared<DFHalfDist>(df->get_imag()->df(), df->get_imag()->nocc());
  const int n = df->get_real()->block().size();
  for (int i=0; i!=n; ++i) {
    dfri->add_block(df->get_real()->block(i)->copy());
    dfri->block(i)->ax_plus_y(1.0, df->get_imag()->block(i));
  }
  auto ricoeff = make_shared<Matrix>(*rcoeff[index] + *icoeff[index]);

  dffull_[0] = df->get_real()->compute_second_transform(rcoeff[index]);
  auto tmp = df->get_imag()->compute_second_transform(icoeff[index]);
  dffull_[1] = dfri->compute_second_transform(ricoeff);
  dffull_[1]->ax_plus_y(-1.0, dffull_[0]);
  dffull_[1]->ax_plus_y(-1.0, tmp);
  dffull_[0]->ax_plus_y(-1.0, tmp);
}


RelDFFull::RelDFFull(array<shared_ptr<DFFullDist>,2> a, pair<int,int> cartesian, vector<shared_ptr<const SpinorInfo>> basis) : RelDFBase(cartesian) {
  basis_ = basis;
  dffull_ = a;
}


RelDFFull::RelDFFull(const RelDFFull& o) : RelDFBase(o.cartesian_) {
  basis_ = o.basis_;
  dffull_[0] = o.dffull_[0]->copy();
  dffull_[1] = o.dffull_[1]->copy();
}


shared_ptr<RelDFFull> RelDFFull::apply_J() const {
  array<shared_ptr<DFFullDist>,2> a{{dffull_[0]->apply_J(), dffull_[1]->apply_J()}};
  return make_shared<RelDFFull>(a, cartesian_, basis_);
}


shared_ptr<RelDFFull> RelDFFull::apply_JJ() const {
  array<shared_ptr<DFFullDist>,2> a{{dffull_[0]->apply_JJ(), dffull_[1]->apply_JJ()}};
  return make_shared<RelDFFull>(a, cartesian_, basis_);
}


shared_ptr<RelDFFull> RelDFFull::clone() const {
  array<shared_ptr<DFFullDist>,2> a{{dffull_[0]->clone(), dffull_[1]->clone()}};
  return make_shared<RelDFFull>(a, cartesian_, basis_);
}


shared_ptr<RelDFFull> RelDFFull::swap() const {
  array<shared_ptr<DFFullDist>,2> a{{dffull_[0]->swap(), dffull_[1]->swap()}};
  a[1]->scale(-1.0);
  return make_shared<RelDFFull>(a, cartesian_, basis_);
}


void RelDFFull::ax_plus_y(const complex<double>& a, const RelDFFull& o) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dffull_[0]->ax_plus_y(fac, o.dffull_[0]);
    dffull_[1]->ax_plus_y(fac, o.dffull_[1]);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dffull_[0]->ax_plus_y(-fac, o.dffull_[1]);
    dffull_[1]->ax_plus_y( fac, o.dffull_[0]);
  } else {
    const double rfac = real(a);
    dffull_[0]->ax_plus_y(rfac, o.dffull_[0]);
    dffull_[1]->ax_plus_y(rfac, o.dffull_[1]);
    const double ifac = imag(a);
    dffull_[0]->ax_plus_y(-ifac, o.dffull_[1]);
    dffull_[1]->ax_plus_y( ifac, o.dffull_[0]);
  }
}


void RelDFFull::scale(complex<double> a) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dffull_[0]->scale(fac);
    dffull_[1]->scale(fac);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dffull_[0]->scale( fac);
    dffull_[1]->scale(-fac);
    std::swap(dffull_[0], dffull_[1]);
  } else {
    throw logic_error("should not happen..");
  }
}


shared_ptr<btas::Tensor3<complex<double>>>
  RelDFFull::get_block(const int ist, const int ii, const int jst, const int jj, const int kst, const int kk) const {

  auto out = make_shared<btas::Tensor3<complex<double>>>(ii, jj, kk);

  shared_ptr<btas::Tensor3<double>> real = get_real()->get_block(ist, ii, jst, jj, kst, kk);
  shared_ptr<btas::Tensor3<double>> imag = get_imag()->get_block(ist, ii, jst, jj, kst, kk);
  auto r = real->data();
  auto i = imag->data();
  for (auto& o : *out)
    o = complex<double>(*r++, *i++);

  return out;
}


list<shared_ptr<RelDFHalfB>> RelDFFull::back_transform(array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff) const {
  list<shared_ptr<RelDFHalfB>> out;
  assert(basis_.size() == 1);
  const int alpha = basis_[0]->alpha_comp();

  for (int i = 0; i != 4; ++i) {
    // Note that icoeff should be scaled by -1.0 !!

    shared_ptr<DFHalfDist> real = dffull_[0]->back_transform(rcoeff[i]);
    real->ax_plus_y( 1.0, dffull_[1]->back_transform(icoeff[i]));

    shared_ptr<DFHalfDist> imag = dffull_[1]->back_transform(rcoeff[i]);
    imag->ax_plus_y(-1.0, dffull_[0]->back_transform(icoeff[i]));

    out.push_back(make_shared<RelDFHalfB>(array<shared_ptr<DFHalfDist>,2>{{real, imag}}, i, alpha));
  }
  return out;
}


shared_ptr<ZMatrix> RelDFFull::form_4index_1fixed(shared_ptr<const RelDFFull> a, const double fac, const int i) const {
#ifndef NDEBUG
  const size_t size = dffull_[0]->nocc1() * dffull_[0]->nocc2() * a->dffull_[0]->nocc1();
  if (size != dffull_[1]->nocc1() * dffull_[1]->nocc2() * a->dffull_[1]->nocc1() ||
      size != dffull_[0]->nocc1() * dffull_[0]->nocc2() * a->dffull_[1]->nocc1())
    throw logic_error("illegal call of RelDFFull::form_4index_1fixed");
#endif

  shared_ptr<Matrix> real = dffull_[0]->form_4index_1fixed(a->dffull_[0], fac, i);
  *real += *dffull_[1]->form_4index_1fixed(a->dffull_[1], -fac, i);

  shared_ptr<Matrix> imag = dffull_[0]->form_4index_1fixed(a->dffull_[1], fac, i);
  *imag += *dffull_[1]->form_4index_1fixed(a->dffull_[0], fac, i);

  return make_shared<ZMatrix>(*real, *imag);
}


shared_ptr<ZMatrix> RelDFFull::form_4index(shared_ptr<const RelDFFull> a, const double fac) const {
  shared_ptr<Matrix> real = dffull_[0]->form_4index(a->dffull_[0], fac);
  *real += *dffull_[1]->form_4index(a->dffull_[1], -fac);

  shared_ptr<Matrix> imag = dffull_[0]->form_4index(a->dffull_[1], fac);
  *imag += *dffull_[1]->form_4index(a->dffull_[0], fac);

  return make_shared<ZMatrix>(*real, *imag);
}


shared_ptr<ZMatrix> RelDFFull::form_2index(shared_ptr<const RelDFFull> a, const double fac, const bool conjugate_left) const {
  const double ifac = conjugate_left ? -1.0 : 1.0;

  shared_ptr<Matrix> real = dffull_[0]->form_2index(a->dffull_[0], fac);
  *real -= *dffull_[1]->form_2index(a->dffull_[1], fac*ifac);

  shared_ptr<Matrix> imag = dffull_[0]->form_2index(a->dffull_[1], fac);
  *imag += *dffull_[1]->form_2index(a->dffull_[0], fac*ifac);

  return make_shared<ZMatrix>(*real, *imag);
}


shared_ptr<RelDFFull> RelDFFull::apply_2rdm(shared_ptr<const ZMatrix> rdm2) const {
  const int nocc = nocc1();
  assert(nocc == nocc2() && nocc*nocc == rdm2->ndim() && rdm2->ndim() == rdm2->mdim());
  shared_ptr<Matrix> rrdm = rdm2->get_real_part();
  shared_ptr<Matrix> irdm = rdm2->get_imag_part();
  const btas::CRange<4> range(nocc, nocc, nocc, nocc);
  static_pointer_cast<btas::Tensor2<double>>(rrdm)->resize(range);
  static_pointer_cast<btas::Tensor2<double>>(irdm)->resize(range);

  shared_ptr<DFFullDist> r  =  dffull_[0]->apply_2rdm(*rrdm);
  r->ax_plus_y(-1.0, dffull_[1]->apply_2rdm(*irdm));

  shared_ptr<DFFullDist> i  =  dffull_[1]->apply_2rdm(*rrdm);
  i->ax_plus_y( 1.0, dffull_[0]->apply_2rdm(*irdm));
  return make_shared<RelDFFull>(array<shared_ptr<DFFullDist>,2>{{r, i}}, cartesian_, basis_);
}


shared_ptr<RelDFFull> RelDFFull::apply_2rdm(shared_ptr<const ZRDM<2>> rdm2) const {
  shared_ptr<const RDM<2>> rrdm = rdm2->get_real_part();
  shared_ptr<const RDM<2>> irdm = rdm2->get_imag_part();

  shared_ptr<DFFullDist> r  =  dffull_[0]->apply_2rdm(*rrdm);
  r->ax_plus_y(-1.0, dffull_[1]->apply_2rdm(*irdm));

  shared_ptr<DFFullDist> i  =  dffull_[1]->apply_2rdm(*rrdm);
  i->ax_plus_y( 1.0, dffull_[0]->apply_2rdm(*irdm));
  return make_shared<RelDFFull>(array<shared_ptr<DFFullDist>,2>{{r, i}}, cartesian_, basis_);
}


void ListRelDFFull::ax_plus_y(const complex<double>& a, shared_ptr<const ListRelDFFull> o) {
  assert(data_.size() == o->data_.size());
  auto oi = o->data_.begin();
  for (auto& i : data_)
    i->ax_plus_y(a, *oi++);
}


shared_ptr<ListRelDFFull> ListRelDFFull::copy() const {
  auto out = make_shared<ListRelDFFull>();
  for (auto& i : data_)
    out->push_back(i->copy());
  return out;
}


shared_ptr<ListRelDFFull> ListRelDFFull::clone() const {
  auto out = make_shared<ListRelDFFull>();
  for (auto& i : data_)
    out->push_back(i->clone());
  return out;
}


shared_ptr<ListRelDFFull> ListRelDFFull::swap() const {
  auto out = make_shared<ListRelDFFull>();
  for (auto& i : data_)
    out->push_back(i->swap());
  return out;
}


shared_ptr<ListRelDFFull> ListRelDFFull::apply_2rdm(shared_ptr<const ZRDM<2>> rdm2) const {
  assert(data_.front()->nocc1() == data_.front()->nocc2() && data_.front()->nocc1() == rdm2->norb());
  auto out = make_shared<ListRelDFFull>();
  for (auto& i : data_)
    out->push_back(i->apply_2rdm(rdm2));
  return out;
}


shared_ptr<ListRelDFFull> ListRelDFFull::apply_2rdm(shared_ptr<const ZMatrix> rdm2) const {
  const int norb = data_.front()->nocc1();
  assert(norb == data_.front()->nocc2() && norb*norb == rdm2->ndim() && norb*norb == rdm2->mdim());
  auto rdm2c = make_shared<ZRDM<2>>(norb);
  copy_n(rdm2->data(), rdm2->size(), rdm2c->data());
  return apply_2rdm(rdm2c);
}


shared_ptr<ZMatrix> ListRelDFFull::form_4index(shared_ptr<const ListRelDFFull> o, const double fac) const {
  shared_ptr<ZMatrix> out;
  for (auto& ii : data_)
    for (auto& jj : o->data_)
      if (ii->alpha_matches(jj)) {
        if (out) {
          *out += *ii->form_4index(jj, fac);
        } else {
          out = ii->form_4index(jj, fac);
        }
      }
  assert(out);
  return out;
}


shared_ptr<ZMatrix> ListRelDFFull::form_4index_1fixed(shared_ptr<const ListRelDFFull> o, const double fac, const int i) const {
  shared_ptr<ZMatrix> out;
  for (auto& ii : data_)
    for (auto& jj : o->data_)
      if (ii->alpha_matches(jj)) {
        if (out) {
          *out += *ii->form_4index_1fixed(jj, fac, i);
        } else {
          out = ii->form_4index_1fixed(jj, fac, i);
        }
      }
  assert(out);
  return out;
}


shared_ptr<ZMatrix> ListRelDFFull::form_2index(shared_ptr<const ListRelDFFull> o, const double fac, const bool conjugate_left) const {
  shared_ptr<ZMatrix> out;
  for (auto& ii : data_)
    for (auto& jj : o->data_)
      if (ii->alpha_matches(jj)) {
        if (out) {
          *out += *ii->form_2index(jj, fac, conjugate_left);
        } else {
          out = ii->form_2index(jj, fac, conjugate_left);
        }
      }
  assert(out);
  return out;
}

