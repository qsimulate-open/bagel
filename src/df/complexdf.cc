//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexdf.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/df/complexdf.h>
#include <src/integral/rys/eribatch.h>

using namespace std;
using namespace bagel;


// Default constructor
ComplexDFDist::ComplexDFDist(const int nbas, const int naux, const array<shared_ptr<DFBlock>,2> blocks, shared_ptr<const ParallelDF> df, shared_ptr<Matrix> data2)
  : DFDist (nbas, naux, nullptr, df, data2), ComplexDF_base() {
  assert((blocks[0] && blocks[1]) || (!blocks[0] && !blocks[1]));
  if (blocks[0]) {
    block_.push_back(blocks[0]);
    block_.push_back(blocks[1]);
  }
  assign_complex_blocks(*this);
}


vector<shared_ptr<const DFDist>> ComplexDFDist::split_blocks() const {
  assert(nindex1_ == nindex2_);
  assert(block_.size() % 2 == 0);
  vector<shared_ptr<const DFDist>> out;
  for (int i = 0; i != block_.size()/2; ++i) {
    array<shared_ptr<DFBlock>,2> in = {{ block_[2*i], block_[2*i+1] }};
    out.push_back(make_shared<const ComplexDFDist>(nindex1_, naux_, in, df_, data2_));
  }
  return out;
}


array<shared_ptr<const DFDist>,2> ComplexDFDist::split_real_imag() const {
  assert(nindex1_ == nindex2_);
  assert(block_.size() % 2 == 0);
  auto outr = make_shared<DFDist>(nindex1_, naux_, nullptr, df_, data2_);
  auto outi = make_shared<DFDist>(nindex1_, naux_, nullptr, df_, data2_);
  for (int i = 0; i != block_.size()/2; ++i) {
    outr->add_block(block_[2*i]);
    outi->add_block(block_[2*i+1]);
  }
  array<shared_ptr<const DFDist>,2> out = {{outr, outi}};
  return out;
}


// Note that we are transforming the bra index, so we need the complex conjugate
shared_ptr<ComplexDFHalfDist> ComplexDFDist::complex_compute_half_transform(const ZMatView c) const {

  // TODO Avoid unnecessary copying here
  auto cmat = make_shared<const ZMatrix>(c);
  shared_ptr<const Matrix> cr = cmat->get_real_part();
  shared_ptr<const Matrix> ci = cmat->get_imag_part();
  auto cri = make_shared<Matrix>(*cr - *ci);

  const int nocc = c.extent(1);
  auto out = make_shared<ComplexDFHalfDist>(df_ ? df_ : shared_from_this(), nocc);
  const int n = block_.size() / 2;
  assert(2*n == block_.size());

  for (int i = 0; i != n; ++i) {
    shared_ptr<DFBlock> blockri = block_[2*i]->copy();
    *blockri += *block_[2*i+1];

    auto rpart = block_[2*i]->transform_second(*cr);
    auto tmp = block_[2*i+1]->transform_second(*ci);
    auto ipart = blockri->transform_second(*cri);

    *ipart -= *rpart;
    *ipart += *tmp;
    *rpart += *tmp;

    out->add_block(rpart);
    out->add_block(ipart);
  }

  out->assign_complex_blocks(*out);
  return out;
}


shared_ptr<ComplexDFHalfDist> ComplexDFDist::complex_compute_half_transform_swap(const ZMatView c) const {
  throw runtime_error("ComplexDFDist::complex_compute_half_transform_swap has not been verified - use caution.");

  // TODO Avoid unnecessary copying here
  auto cmat = make_shared<const ZMatrix> (c);
  const shared_ptr<Matrix> cr = cmat->get_real_part();
  const shared_ptr<Matrix> ci = cmat->get_imag_part();
  auto cri = make_shared<Matrix>(*cr + *ci);

  const int nocc = c.extent(1);
  auto out = make_shared<ComplexDFHalfDist>(df_ ? df_ : shared_from_this(), nocc);
  const int n = block_.size() / 2;
  assert(2*n == block_.size());

  for (int i = 0; i != n; ++i) {
    shared_ptr<DFBlock> blockri = block_[2*i]->copy();
    *blockri += *block_[2*i+1];

    auto rpart = block_[2*i]->transform_third(*cr);
    auto tmp = block_[2*i+1]->transform_third(*ci);
    auto ipart = blockri->transform_third(*cri);

    *ipart -= *rpart;
    *ipart -= *tmp;
    *rpart -= *tmp;

    out->add_block(rpart);
    out->add_block(ipart);
  }
  out->assign_complex_blocks(*out);
  return out;
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop_from_cd(shared_ptr<const ZVectorB> tmp0) const {
  if (block_.size() != 2) throw logic_error("compute_Jop so far assumes block_.size() == 2 for complex integrals");
  const shared_ptr<VectorB> dr = tmp0->get_real_part();
  const shared_ptr<VectorB> di = tmp0->get_imag_part();
  shared_ptr<Matrix> outr, outi;

  {
    auto dri = make_shared<VectorB>(*dr + *di);
    shared_ptr<DFBlock> blockri = block_[0]->copy();
    *blockri += *block_[1];

    outr = block_[0]->form_mat(dr->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize()));
    outi = blockri->form_mat(dri->slice(blockri->astart(), blockri->astart()+blockri->asize()));
    shared_ptr<Matrix> tmp = block_[1]->form_mat(di->slice(block_[1]->astart(), block_[1]->astart()+block_[1]->asize()));

    *outi -= *outr;
    *outi -= *tmp;
    *outr -= *tmp;
  }

  auto out = make_shared<ZMatrix>(*outr, *outi);

  // all reduce
  if (!serial_)
    out->allreduce();

  return out;
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop(const shared_ptr<const ZMatrix> den) const {
  auto tmp = dynamic_pointer_cast<const ComplexDFDist>(shared_from_this());
  assert(tmp);
  return complex_compute_Jop(tmp, den);
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop(const shared_ptr<const ComplexDF_base> o, const shared_ptr<const ZMatrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const ZVectorB> tmp0 = o->complex_compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return complex_compute_Jop_from_cd(tmp0);
}


shared_ptr<ComplexDFHalfDist> ComplexDFHalfDist::complex_apply_J(const shared_ptr<const Matrix> d) const {
  auto tmp = apply_J(d);
  auto out = make_shared<ComplexDFHalfDist>(tmp, nindex1_);
  out->block_ = tmp->block();
  out->assign_complex_blocks(*out);
  return out;
}


shared_ptr<ZMatrix> ComplexDFHalfDist::complex_form_2index(shared_ptr<const ComplexDFHalfDist> o, const double a, const bool swap) const {
  if (block_.size() != 2 || o->block_.size() != 2) throw logic_error("so far assumes block_.size() == 2 for complex integrals");

  shared_ptr<ZMatrix> out;

  if (!swap) {
    shared_ptr<DFBlock> blockri = block_[0]->copy();
    *blockri -= *block_[1];
    shared_ptr<DFBlock> oblockri = o->block_[0]->copy();
    *oblockri += *o->block_[1];

    shared_ptr<Matrix> rout = block_[0]->form_2index(o->block_[0], a);
    shared_ptr<Matrix> tmp = block_[1]->form_2index(o->block_[1], a);
    shared_ptr<Matrix> iout = blockri->form_2index(oblockri, a);

    *iout -= *rout;
    *iout += *tmp;
    *rout += *tmp;
    out = make_shared<ZMatrix>(*rout, *iout);
  } else {
    throw runtime_error("Please verify carefully that ComplexDFHalfDist::form_2index(...) is correct when swap = true");
#if 0
    shared_ptr<DFBlock> oblockri = o->block_[0]->copy();
    *oblockri -= *o->block_[1];
    shared_ptr<DFBlock> blockri = block_[0]->copy();
    *blockri += *block_[1];

    shared_ptr<Matrix> rout = o->block_[0]->form_2index(block_[0], a);
    shared_ptr<Matrix> tmp = o->block_[1]->form_2index(block_[1], a);
    shared_ptr<Matrix> iout = oblockri->form_2index(blockri, a);

    *iout -= *rout;
    *iout += *tmp;
    *rout += *tmp;

    out = make_shared<ZMatrix>(*rout, *iout);
#endif
  }

  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<ComplexDFFullDist> ComplexDFHalfDist::complex_compute_second_transform(const ZMatView c) const {
  // TODO Avoid unnecessary copying here
  auto cmat = make_shared<const ZMatrix>(c);
  shared_ptr<const Matrix> cr = cmat->get_real_part();
  shared_ptr<const Matrix> ci = cmat->get_imag_part();

  auto out = make_shared<ComplexDFFullDist>(df_ ? df_ : shared_from_this(), nocc(), c.extent(1));
  const int n = block_.size() / 2;
  assert(2*n == block_.size());

  for (int i = 0; i != n; ++i) {
    auto rpart =  block_[2*i]->transform_third(*cr);
    *rpart    -= *block_[2*i+1]->transform_third(*ci);
    out->add_block(rpart);

    auto ipart =  block_[2*i+1]->transform_third(*cr);
    *ipart    += *block_[2*i]->transform_third(*ci);
    out->add_block(ipart);
  }

  out->assign_complex_blocks(*out);
  return out;
}


shared_ptr<ZMatrix> ComplexDFFullDist::complex_form_4index(shared_ptr<const ComplexDFFullDist> o, const double a) const {
  if (block_.size() != 2 || o->block_.size() != 2)
    throw logic_error("so far ComplexDFFullDist::complex_form_4index assumes block_.size() == 1");

  auto rpart =  block_[0]->form_4index(o->block_[0], a);
  *rpart    -= *block_[1]->form_4index(o->block_[1], a);

  auto ipart =  block_[0]->form_4index(o->block_[1], a);
  *ipart    += *block_[1]->form_4index(o->block_[0], a);

  return make_shared<ZMatrix>(*rpart, *ipart);
}
