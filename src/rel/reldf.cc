//
// BAGEL - Parallel electron correlation program.
// Filename: reldf.cc
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


#include <src/rel/reldf.h>
#include <src/df/complexdf.h>

using namespace std;
using namespace bagel;

RelDF::RelDF(shared_ptr<const DFDist> df, pair<int, int> cartesian, const std::vector<int> alpha) : RelDFBase(cartesian), alpha_(alpha), swap_(false) {
  if (!dynamic_pointer_cast<const ComplexDFDist>(df)) {
    dfdata_[0] = df;
    dfdata_[1] = nullptr;
  } else {
    auto cdf = dynamic_pointer_cast<const ComplexDFDist>(df);
    array<shared_ptr<const DFDist>,2> tmp = cdf->split_real_imag();
    dfdata_[0] = tmp[0];
    dfdata_[1] = tmp[1];
  }
  set_basis();
}

RelDF::RelDF(const RelDF& o, bool coo) : RelDFBase(o), alpha_(o.alpha_), dfdata_(o.get_data()), swap_(o.swap_) {
  set_basis();
  if (coo) {
    swap_ ^= true;
    vector<std::shared_ptr<const SpinorInfo>> newbas;
    for (auto& i : basis_)
      newbas.push_back(i->swap());
    basis_ = newbas;
    std::swap(cartesian_.first, cartesian_.second);
  }
}


//swap cartesian
shared_ptr<const RelDF> RelDF::swap() const {
  return make_shared<const RelDF>(*this, true);
}


vector<shared_ptr<RelDFHalf>> RelDF::compute_half_transform(array<shared_ptr<const Matrix>,4> rc, array<shared_ptr<const Matrix>,4> ic) const {
  vector<shared_ptr<RelDFHalf>> out;

  // first make a subset
  vector<vector<shared_ptr<const SpinorInfo>>> subsets;
  for (int i = 0; i != 4; ++i) {
    vector<shared_ptr<const SpinorInfo>> tmp;
    for (auto& j : basis())
      if (j->basis(0) == i)
        tmp.push_back(j);
    if (!tmp.empty())
      subsets.push_back(tmp);
  }

  // transform
  for (auto& i : subsets)
    out.push_back(make_shared<RelDFHalf>(shared_from_this(), i, rc, ic));
  return out;
}


vector<shared_ptr<ZMatrix>> RelDF::compute_Jop(list<shared_ptr<const CDMatrix>>& cd) const {

  vector<shared_ptr<ZVectorB>> sum;
  for (auto& b : basis_) {
    sum.push_back(cd.front()->clone());
    for (auto& i : cd) {
      if(b->alpha_comp() == i->alpha_comp())
        *sum.back() += *i;
    }
  }

  vector<shared_ptr<ZMatrix>> out;

  // real
  if (!get_imag()) {
    for (auto& i : sum) {
      shared_ptr<const Matrix> rdat = get_real()->compute_Jop_from_cd(i->get_real_part());
      shared_ptr<const Matrix> idat = get_real()->compute_Jop_from_cd(i->get_imag_part());
      out.push_back(make_shared<ZMatrix>(*rdat, *idat));
    }

  // complex (GIAO)
  } else {
    // get_real() + get_imag() needed for 3-multiplication algorithm
    auto dfri = make_shared<DFDist>(get_real());
    const int n = get_real()->block().size();
    for (int i=0; i!=n; ++i) {
      dfri->add_block(get_real()->block(i)->copy());
      dfri->block(i)->ax_plus_y(1.0, get_imag()->block(i));
    }
    for (auto& i : sum) {
      shared_ptr<const VectorB> cdr = i->get_real_part();
      shared_ptr<const VectorB> cdi = i->get_imag_part();
      auto cdri = make_shared<const VectorB> (*cdr + *cdi);

      shared_ptr<Matrix> rdat = get_real()->compute_Jop_from_cd(cdr);
      shared_ptr<Matrix> tmpdat = get_imag()->compute_Jop_from_cd(cdi);
      shared_ptr<Matrix> idat = dfri->compute_Jop_from_cd(cdri);
      *idat -= *rdat;
      *idat -= *tmpdat;
      *rdat -= *tmpdat;
      out.push_back(make_shared<ZMatrix>(*rdat, *idat));
    }
  }

  return out;
}


