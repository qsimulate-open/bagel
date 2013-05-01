//
// BAGEL - Parallel electron correlation program.
// Filename: relfci.cc
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
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

#include <string>
#include <vector>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <src/rel/relfci.h>
#include <src/fci/space.h>
#include <src/rysint/eribatch.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <src/util/davidson.h>

using namespace std;
using namespace bagel;

RelFCI::RelFCI(const boost::property_tree::ptree& idata, const shared_ptr<const Geometry> geom,
             const shared_ptr<const RelReference> re) : relref_(re) {
  gaunt_ = idata.get<bool>("gaunt", true);
  breit_ = idata.get<bool>("breit", gaunt_);
  geom_ = geom->relativistic(gaunt_);
  common_init(idata);
}


void RelFCI::common_init(const boost::property_tree::ptree& idata) {
  cout << "  *** Relativistic Full CI ***" << endl << endl;

  // reading input keywords
  max_iter_ = idata.get<int>("maxiter", 100);
  max_iter_ = idata.get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata.get<int>("diis_start", 1);
  thresh_scf_ = idata.get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata.get<double>("thresh_scf", thresh_scf_);
  ncharge_ = idata.get<int>("charge", 0);
  nele_ = geom_->nele()-ncharge_;
  nneg_ = geom_->nbasis()*2;

  if (breit_ && !gaunt_) throw runtime_error("Breit cannot be turned on if Gaunt is off");
}


shared_ptr<const ZMatrix> RelFCI::time_reversal_operator() {

  const int n = geom_->nbasis();
  std::complex<double> one  (1.0, 0.0);
  std::complex<double> coeffi  (0.0, 1.0);
  
  auto kramers = make_shared<ZMatrix>(4*n, 4*n);

  auto unit = make_shared<ZMatrix>(n, n);
  unit->unit();

  kramers->add_block(-1, 0, n, n, n, unit);
  kramers->add_block( 1, n, 0, n, n, unit);
  kramers->add_block(-1, 2*n, 3*n, n, n, unit);
  kramers->add_block( 1, 3*n, 2*n, n, n, unit);

  auto overlap = make_shared<RelOverlap>(geom_, false);

  return make_shared<const ZMatrix>(*kramers * *overlap * coeffi);
}


void RelFCI::print_eig(const unique_ptr<double[]>& eig, const int n) {
  for (int i = 0; i != n; ++i) cout << setprecision(10) << setw(15) << eig[i] << "    " <<  i << endl;
}


void RelFCI::compute() {
  cout << "Debug print out in RelFCI::compute()... " << endl;
  unique_ptr<double[]> eig(new double[relref_->coeff()->ndim()]);
  unique_ptr<double[]> eig2(new double[relref_->coeff()->ndim()]);

  coeff_ = relref_->coeff();
  shared_ptr<const ZMatrix> coeff = coeff_;
  nocc_ = relref_->nocc();
  nvirt_ = relref_->nvirt();
  shared_ptr<const ZMatrix> coeff_occ = coeff->slice(0, nocc_);
  shared_ptr<const ZMatrix> coeff_virt = coeff->slice(nocc_, nocc_+nvirt_);
  shared_ptr<const ZMatrix> coeff_occ_conjg = coeff_occ->get_conjg();
  shared_ptr<const ZMatrix> coeff_virt_conjg = coeff_virt->get_conjg();

  shared_ptr<ZMatrix> time_reversal(make_shared<ZMatrix>(*time_reversal_operator()));
#if 1
  ZMatrix tmp1(*coeff_occ % *time_reversal * *coeff_occ_conjg);
  ZMatrix tmp2(*coeff_virt % *time_reversal * *coeff_virt_conjg);
#else
  ZMatrix tmp1(*coeff_occ % *time_reversal * *coeff_occ);
  ZMatrix tmp2(*coeff_virt % *time_reversal * *coeff_virt);
#endif
tmp1.print("T");
#if 0
  tmp1.diagonalize_skew(eig.get());
  tmp2.diagonalize_skew(eig2.get());
#else
  tmp1.diagonalize(eig.get());
  tmp2.diagonalize(eig2.get());
#endif
  cout << " OCCUPIED " << endl;
  print_eig(eig, tmp1.mdim());
  cout << " VIRTUAL " << endl;
  print_eig(eig2, tmp2.mdim());

}

