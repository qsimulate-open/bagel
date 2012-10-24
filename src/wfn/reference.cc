//
// BAGEL - Parallel electron correlation program.
// Filename: reference.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/util/f77.h>
#include <src/wfn/reference.h>
#include <src/fci/knowles.h>
#include <src/osint/overlapbatch.h>
#include <src/util/mixedbasis.h>

using namespace std;
using namespace bagel;

Reference::Reference(shared_ptr<const Geometry> g, shared_ptr<const Coeff> c,
                     const int _nclosed, const int _nact, const int _nvirt,
                     const double en,
                     const vector<shared_ptr<RDM<1> > >& _rdm1, const vector<shared_ptr<RDM<2> > >& _rdm2,
                     shared_ptr<const RDM<1> > _rdm1_av, shared_ptr<const RDM<2> > _rdm2_av)
 : geom_(g), coeff_(c), energy_(en), hcore_(new Hcore(geom_)), nclosed_(_nclosed), nact_(_nact), nvirt_(_nvirt), nstate_(1), rdm1_(_rdm1), rdm2_(_rdm2),
   rdm1_av_(_rdm1_av), rdm2_av_(_rdm2_av) {

  if (nact_ && rdm1_.empty())
    throw logic_error("If nact != 0, Reference::Reference wants to have RDMs.");

}


shared_ptr<Matrix1e> Reference::rdm1_mat(shared_ptr<const RDM<1> > active) const {
  if (nact_)
    return active->rdm1_mat(geom_, nclosed_);
  else {
    shared_ptr<Matrix1e> out(new Matrix1e(geom_, nocc(), nocc()));
    for (int i = 0; i != nclosed_; ++i) out->element(i,i) = 2.0;
    return out;
  }
}


shared_ptr<Dvec> Reference::civectors() const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci(new KnowlesHandy(multimap<string, string>(), shared_from_this(), nclosed_, nact_, nstate_));
  fci->compute();
  return fci->civectors();
}


// TODO this is a very bad implementation, since it recomputes FCI; should be replaced in somewhere.
shared_ptr<RDM<3> > Reference::compute_rdm3(const int i) const {
  // Default to HarrisonZarrabian method
  shared_ptr<FCI> fci(new KnowlesHandy(multimap<string, string>(), shared_from_this(), nclosed_, nact_, nstate_));
  fci->compute();
  fci->compute_rdm12();
  return get<0>(fci->compute_rdm34(i));
}


shared_ptr<const Reference> Reference::project_coeff(shared_ptr<const Geometry> geomin) const {
  shared_ptr<const Reference> out;

  if (*geom_ == *geomin) {
    out = shared_from_this();
  } else {
    // in this case we first form overlap matrices
    shared_ptr<Overlap> snew(new Overlap(geomin));
    shared_ptr<Overlap> snew2(new Overlap(*snew));
    shared_ptr<MixedBasis<OverlapBatch> > mixed(new MixedBasis<OverlapBatch>(geom_, geomin));
    const int nnew = geomin->nbasis();
    const int nold = geom_->nbasis();

    unique_ptr<int[]> ipiv(new int[nnew+1]);
    dgesv_(nnew, nold, snew->data(), nnew, ipiv.get(), mixed->data(), nnew, ipiv[nnew]);
    if (ipiv[nnew]) throw runtime_error("DGESV failed in Reference::project_coeff");

    shared_ptr<Coeff> c(new Coeff(geomin));
    dgemm_("N", "N", nnew, nold, nold, 1.0, mixed->data(), nnew, coeff_->data(), nold, 0.0, c->data(), nnew);

#if 1
    unique_ptr<double[]> diag(new double[nnew]);
    Matrix1e m = *c % *snew2 * *c;
    for (int i = nold; i < nnew; ++i) m.element(i,i) = 1.0;
    m.diagonalize(diag.get());
    for (int i = 0; i != nnew; ++i) dscal_(nnew, 1.0/sqrt(sqrt(diag[i])), m.data()+i*nnew, 1);
    *c *= (m ^ m);
#endif

    out = shared_ptr<const Reference>(new Reference(geomin, c, nclosed_, nact_, geomin->nbasis()-nclosed_-nact_, energy_));
  }

  return out;
}
