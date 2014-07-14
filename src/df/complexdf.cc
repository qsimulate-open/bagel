//
// BAGEL - Parallel electron correlation program.
// Filename: complexdf.cc
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


#include <src/df/complexdf.h>
#include <src/integral/rys/eribatch.h>

using namespace std;
using namespace bagel;


// Default constructor
ComplexDFDist::ComplexDFDist(const int nbas, const int naux, const array<shared_ptr<DFBlock>,2> blocks, shared_ptr<const ParallelDF> df, shared_ptr<Matrix> data2)
  : DFDist (nbas, naux, nullptr, df, data2) {
  assert((blocks[0] && blocks[1]) || (!blocks[0] && !blocks[1]));
  if (blocks[0]) {
    block_.push_back(blocks[0]);
    block_.push_back(blocks[1]);
  }
}


vector<shared_ptr<const DFDist>> ComplexDFDist::split_complex_blocks() const {
  assert(nindex1_ == nindex2_);
  assert(block_.size() % 2 == 0);
  vector<shared_ptr<const DFDist>> out;
  for (int i=0; i<block_.size(); i+=2) {
    array<shared_ptr<DFBlock>,2> in = {{ block_[i], block_[i+1] }};
    out.push_back(make_shared<const ComplexDFDist>(nindex1_, naux_, in, df_, data2_));
  }
  return out;
}


array<shared_ptr<const DFDist>,2> ComplexDFDist::split_real_imag() const {
  assert(nindex1_ == nindex2_);
  assert(block_.size() % 2 == 0);
  auto outr = make_shared<DFDist>(nindex1_, naux_, nullptr, df_, data2_);
  auto outi = make_shared<DFDist>(nindex1_, naux_, nullptr, df_, data2_);
  for (int i=0; i<block_.size(); i+=2) {
    outr->add_block(block_[i]);
    outi->add_block(block_[i+1]);
  }
  array<shared_ptr<const DFDist>,2> out = {{outr, outi}};
  return out;
}


// Note that we are transforming the bra index, so we need the complex conjugate
shared_ptr<ComplexDFHalfDist> ComplexDFDist::complex_compute_half_transform(const ZMatView c) const {

  // TODO Avoid unnecessary copying here
  auto cmat = make_shared<const ZMatrix> (c);
  const shared_ptr<Matrix> cr = cmat->get_real_part();
  const shared_ptr<Matrix> ci = cmat->get_imag_part();

  const int nocc = c.extent(1);
  auto out = make_shared<ComplexDFHalfDist>(shared_from_this(), nocc);
  const int n = block_.size() / 2;
  assert(2*n == block_.size());

  // TODO Using 4-multiplication
  for (int i=0; i!=n; ++i) {
    auto rpart = block_[i]->transform_second(*cr);
    rpart->ax_plus_y( 1.0, block_[i+1]->transform_second(*ci));
    auto ipart = block_[i+1]->transform_second(*cr);
    ipart->ax_plus_y(-1.0, block_[i]->transform_second(*ci));
    out->add_block(rpart);
    out->add_block(ipart);
  }

  return out;
}


shared_ptr<ComplexDFHalfDist> ComplexDFDist::complex_compute_half_transform_swap(const ZMatView c) const {
  throw runtime_error("ComplexDFDist::compute_half_transform_swap has not been verified - use caution.");

  // TODO Avoid unnecessary copying here
  auto cmat = make_shared<const ZMatrix> (c);
  const shared_ptr<Matrix> cr = cmat->get_real_part();
  const shared_ptr<Matrix> ci = cmat->get_imag_part();

  const int nocc = c.extent(1);
  auto out = make_shared<ComplexDFHalfDist>(shared_from_this(), nocc);
  const int n = block_.size() / 2;
  assert(2*n == block_.size());

  // TODO Using 4-multiplication
  for (int i=0; i!=n; ++i) {
    auto rpart = block_[i]->transform_second(*cr);
    rpart->ax_plus_y(-1.0, block_[i+1]->transform_third(*ci));
    auto ipart = block_[i+1]->transform_second(*cr);
    ipart->ax_plus_y( 1.0, block_[i]->transform_third(*ci));
    out->add_block(rpart);
    out->add_block(ipart);
  }
  return out;
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop_from_cd(shared_ptr<const ZVectorB> tmp0) const {
  if (block_.size() != 2) throw logic_error("compute_Jop so far assumes block_.size() == 2 for complex integrals");
  const shared_ptr<VectorB> dr = tmp0->get_real_part();
  const shared_ptr<VectorB> di = tmp0->get_imag_part();

  // TODO using 4-multiplication
  shared_ptr<Matrix> outr = block_[0]->form_mat(*dr->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize()));
  shared_ptr<Matrix> outi = block_[0]->form_mat(*di->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize()));
  *outr -= *block_[1]->form_mat(*di->slice(block_[1]->astart(), block_[1]->astart()+block_[1]->asize()));
  *outi += *block_[1]->form_mat(*dr->slice(block_[1]->astart(), block_[1]->astart()+block_[1]->asize()));

  // all reduce
  if (!serial_)
    outr->allreduce();

  auto out = make_shared<ZMatrix>(*outr, *outi);
  return out;
}


shared_ptr<ZVectorB> ComplexDFDist::complex_compute_cd(const shared_ptr<const ZMatrix> den, shared_ptr<const Matrix> dat2, const bool onlyonce) const {
  if (block_.size() != 2) throw logic_error("compute_Jop so far assumes block_.size() == 2 for complex integrals");
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = make_shared<ZVectorB>(naux_);
  const shared_ptr<Matrix> dr = den->get_real_part();
  const shared_ptr<Matrix> di = den->get_imag_part();

  // TODO Using 4-multiplication
  // D = (D|rs)*d_rs
  shared_ptr<VectorB> tmpr = block_[0]->form_vec(dr);
  *tmpr -= *block_[1]->form_vec(di);
  shared_ptr<VectorB> tmpi = block_[0]->form_vec(di);
  *tmpi += *block_[1]->form_vec(dr);

  auto tmp = make_shared<ZVectorB>(*tmpr, *tmpi);
  copy_n(tmp->data(), block_[0]->asize(), tmp0->data()+block_[0]->astart());

  // All reduce
  if (!serial_)
    tmp0->allreduce();

  // TODO Inefficient.  Multiply by dat2 before combining real and complex
  auto cdat2 = make_shared<ZMatrix>(*dat2, 1.0);
  *tmp0 = *cdat2 * *tmp0;
  if (!onlyonce)
    *tmp0 = *cdat2 * *tmp0;
  return tmp0;
}


shared_ptr<ZVectorB> ComplexDFHalfDist::complex_compute_cd(const shared_ptr<const ZMatrix> den, shared_ptr<const Matrix> dat2, const bool onlyonce) const {
  if (block_.size() != 2) throw logic_error("compute_Jop so far assumes block_.size() == 2 for complex integrals");
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = make_shared<ZVectorB>(naux_);
  const shared_ptr<Matrix> dr = den->get_real_part();
  const shared_ptr<Matrix> di = den->get_imag_part();

  // TODO Using 4-multiplication
  // D = (D|rs)*d_rs
  shared_ptr<VectorB> tmpr = block_[0]->form_vec(dr);
  *tmpr -= *block_[1]->form_vec(di);
  shared_ptr<VectorB> tmpi = block_[0]->form_vec(di);
  *tmpi += *block_[1]->form_vec(dr);

  auto tmp = make_shared<ZVectorB>(*tmpr, *tmpi);
  copy_n(tmp->data(), block_[0]->asize(), tmp0->data()+block_[0]->astart());

  // All reduce
  if (!serial_)
    tmp0->allreduce();

  // TODO Inefficient.  Multiply by dat2 before combining real and complex
  auto cdat2 = make_shared<ZMatrix>(*dat2, 1.0);
  *tmp0 = *cdat2 * *tmp0;
  if (!onlyonce)
    *tmp0 = *cdat2 * *tmp0;
  return tmp0;
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop(const shared_ptr<const ZMatrix> den) const {
  auto tmp = dynamic_pointer_cast<const ComplexDFDist>(this->shared_from_this());
  assert(tmp);
  return complex_compute_Jop(tmp, den);
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop(const shared_ptr<const ComplexDFDist> o, const shared_ptr<const ZMatrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const ZVectorB> tmp0 = o->complex_compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return complex_compute_Jop_from_cd(tmp0);
}


shared_ptr<ZMatrix> ComplexDFDist::complex_compute_Jop(const shared_ptr<const ComplexDFHalfDist> o, const shared_ptr<const ZMatrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const ZVectorB> tmp0 = o->complex_compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return complex_compute_Jop_from_cd(tmp0);
}


// TODO Essentially the same as DFHalfDist::apply_J, except for the return type
shared_ptr<ComplexDFHalfDist> ComplexDFHalfDist::complex_apply_J(const shared_ptr<const Matrix> d) const {
  auto out = make_shared<ComplexDFHalfDist>(df_, nindex1_);
  for (auto& i : block_)
    out->add_block(i->clone());
#ifdef HAVE_MPI_H
  if (!serial_) {
    Timer mult(3);
    auto work = make_shared<DFDistT>(shared_from_this());
    mult.tick_print("Form DFDistT");
    work = work->apply_J(d);
    mult.tick_print("Application of Inverse");
    work->get_paralleldf(out);
    mult.tick_print("Return DFDist");
  } else {
#else
  {
#endif
    auto j = block_.begin();
    for (auto& i : out->block_) {
      i->zero();
      i->contrib_apply_J(*j, d);
      ++j;
    }
  }
  return out;
}


shared_ptr<ZMatrix> ComplexDFHalfDist::complex_form_2index(shared_ptr<const ComplexDFHalfDist> o, const double a, const bool swap) const {
  if (block_.size() != 2 || o->block_.size() != 2) throw logic_error("so far assumes block_.size() == 2 for complex integrals");

  shared_ptr<ZMatrix> out;

  // TODO using 4-multiplication
  if (!swap) {
    shared_ptr<Matrix> rout = block_[0]->form_2index(o->block_[0], a);
    shared_ptr<Matrix> iout = block_[0]->form_2index(o->block_[1], a);
    *rout += *block_[1]->form_2index(o->block_[1], a);
    *iout -= *block_[1]->form_2index(o->block_[0], a);
    out = make_shared<ZMatrix>(*rout, *iout);
  } else {
    throw runtime_error("Please verify carefully that ComplexDFHalfDist::form_2index(...) is correct when swap = true");
    shared_ptr<Matrix> rout = o->block_[0]->form_2index(block_[0], a);
    shared_ptr<Matrix> iout = o->block_[0]->form_2index(block_[1], a);
    *rout += *o->block_[1]->form_2index(block_[1], a);
    *iout -= *o->block_[1]->form_2index(block_[0], a);
    out = make_shared<ZMatrix>(*rout, *iout);
  }
  if (!serial_)
    out->allreduce();
  return out;
}

