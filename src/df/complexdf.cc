//
// BAGEL - Parallel electron correlation program.
// Filename: complex.cc
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
ComplexDFDist::ComplexDFDist(const int nbas, const int naux, const int nblock, const shared_ptr<DFBlock> block, shared_ptr<const ComplexParallelDF> df, shared_ptr<Matrix> data2)
  : ComplexParallelDF(naux, nbas, nbas, nblock, df, data2) {
  dfdata_[0] = make_shared<DFDist>(nbas, naux, block, (df ? df->get_real() : nullptr), data2);
  dfdata_[1] = make_shared<DFDist>(nbas, naux, block, (df ? df->get_imag() : nullptr), data2);
  if (block)
    throw std::runtime_error("How do we want to handle the construction of ComplexDFDists, when provided with a DFBlock?");
    //block_.push_back(block);
}


ComplexDFDist::ComplexDFDist(const shared_ptr<const ComplexParallelDF> df) : ComplexParallelDF(df->naux(), df->nindex1(), df->nindex2(), df->nblock(), df) { }


// Constructor used by split_blocks()
ComplexDFDist::ComplexDFDist(const std::array<std::shared_ptr<DFBlock>,2> block, const std::shared_ptr<const ComplexParallelDF> df, const shared_ptr<Matrix> d2)
 : ComplexParallelDF(df->naux(), df->nindex1(), df->nindex2(), 1) {
  assert(block[0]->asize() == block[1]->asize());
  assert(block[0]->b1size() == block[1]->b1size());
  assert(block[0]->b2size() == block[1]->b2size());
  assert(block[0]->b1size() == block[0]->b2size());
  auto dfr = make_shared<DFDist> (nindex1_, naux_, block[0], (df ? df->get_real() : nullptr), d2);
  auto dfi = make_shared<DFDist> (nindex1_, naux_, block[1], (df ? df->get_imag() : nullptr), d2);
  dfdata_ = {{ dfr, dfi }};
  data2_ = d2;
}


shared_ptr<const StaticDist> ComplexDFDist::make_table(const size_t astart) {
  vector<size_t> rec(mpi__->size());
  fill(rec.begin(), rec.end(), 0);

  mpi__->allgather(&astart, 1, rec.data(), 1);
  rec.push_back(naux_);

  return make_shared<const StaticDist>(rec);
}


tuple<int, vector<shared_ptr<const Shell>>> ComplexDFDist::get_ashell(const vector<shared_ptr<const Shell>>& all) {
  int out1;
  vector<shared_ptr<const Shell>> out2;
  // TODO without *2, H does not work. Perhaps need to think a bit more
  if (mpi__->size()*2 < all.size()) {
    int start, end;
    StaticDist d(naux_, mpi__->size());
    tie(start, end) = d.range(mpi__->rank());
    int num = 0;
    for (auto iter = all.begin(); iter != all.end(); ++iter) {
      if (num >= start && num < end) {
        if (out2.empty()) out1 = num;
        out2.push_back(*iter);
      }
      num += (*iter)->nbasis();
    }
  } else {
    cout << endl << "   *** Warning *** Since the number of auxiliary shells is too small, we do not parallelize the Fock builder." << endl << endl;
    out1 = 0;
    out2 = all;
    serial_ = true;
  }

  return tie(out1, out2);
}


void ComplexDFDist::compute_2index(const vector<shared_ptr<const Shell>>& ashell, const double throverlap, const bool compute_inverse) {

  auto tmp = dynamic_pointer_cast<DFDist>(dfdata_[0]);
  if (!tmp) throw runtime_error("Illegal call of ComplexDFDist::compute_2index");

  tmp->compute_2index(ashell, throverlap, compute_inverse);
  dfdata_[1]->data2_ = tmp->data2_;
  data2_ = tmp->data2_;

  assert(dfdata_[0]->data2_ == dfdata_[1]->data2_);
}


// Note that we are transforming the bra index, so we need the complex conjugate
shared_ptr<ComplexDFHalfDist> ComplexDFDist::compute_half_transform(const ZMatView c) const {
  auto cmat = make_shared<const ZMatrix> (c);
  const std::shared_ptr<Matrix> cr = cmat->get_real_part();
  const std::shared_ptr<Matrix> ci = cmat->get_imag_part();
  auto rpart = dynamic_pointer_cast<const DFDist>(dfdata_[0]);
  auto ipart = dynamic_pointer_cast<const DFDist>(dfdata_[1]);
  assert(rpart && ipart);
  auto outr = rpart->compute_half_transform(cr);
  auto outi = ipart->compute_half_transform(cr);
  outr->ax_plus_y( 1.0, ipart->compute_half_transform(ci));
  outi->ax_plus_y(-1.0, rpart->compute_half_transform(ci));
  return make_shared<ComplexDFHalfDist> (outr, outi, shared_from_this());
}


shared_ptr<ComplexDFHalfDist> ComplexDFDist::compute_half_transform_swap(const ZMatView c) const {
  throw runtime_error("ComplexDFDist::compute_half_transform_swap has not been verified - use caution.");
  auto cmat = make_shared<const ZMatrix> (c);
  const std::shared_ptr<Matrix> cr = cmat->get_real_part();
  const std::shared_ptr<Matrix> ci = cmat->get_imag_part();
  auto rpart = dynamic_pointer_cast<const DFDist>(dfdata_[0]);
  auto ipart = dynamic_pointer_cast<const DFDist>(dfdata_[1]);
  assert(rpart && ipart);
  shared_ptr<DFHalfDist> outr = rpart->compute_half_transform_swap(cr);
  shared_ptr<DFHalfDist> outi = rpart->compute_half_transform_swap(ci);
  outr->ax_plus_y(-1.0, ipart->compute_half_transform_swap(ci));
  outi->ax_plus_y( 1.0, ipart->compute_half_transform_swap(cr));
  outi->scale(-1.0);
  return make_shared<ComplexDFHalfDist> (outr, outi, shared_from_this());
}


ComplexDFHalfDist::ComplexDFHalfDist(std::shared_ptr<DFHalfDist> rdf, std::shared_ptr<DFHalfDist> idf, std::shared_ptr<const ComplexParallelDF> df)
 : ComplexParallelDF(rdf->naux(), rdf->nocc(), rdf->nbasis(), df->nblock(), df) {
  assert(rdf->naux() == idf->naux());
  assert(rdf->nocc() == idf->nocc());
  assert(rdf->nbasis() == idf->nbasis());
  assert(rdf->block().size() == idf->block().size());
  assert(rdf->block().size() == nblock_);
  dfdata_[0] = rdf;
  dfdata_[1] = idf;
}


shared_ptr<ComplexDFHalfDist> ComplexDFHalfDist::apply_J(const shared_ptr<const Matrix> d) const {
  auto rpart = dynamic_pointer_cast<const DFHalfDist>(dfdata_[0]);
  auto ipart = dynamic_pointer_cast<const DFHalfDist>(dfdata_[1]);
  assert(rpart && ipart);
  auto rout = rpart->apply_J(d);
  auto iout = ipart->apply_J(d);
  auto out = make_shared<ComplexDFHalfDist>(rout, iout, df_);
  return out;
}

