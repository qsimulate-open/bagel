//
// BAGEL - Parallel electron correlation program.
// Filename: soscf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu> 
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

#include <src/scf/soscf.h>

using namespace std;
using namespace bagel;

SOSCF::SOSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : SCF_base(idata, geom, re) {
  cout << "hello world" << endl;
  cout << setprecision(10) << geom->nbasis() << " " << idata->get<double>("thresh", 0.0) << endl;
  cout << setprecision(10) << geom->nbasis() << " " << idata->get<double>("thresha", 0.0) << endl;
}


void SOSCF::compute() {
  cout << "hello world from compute" << endl;

#if 0
    Fock(const std::shared_ptr<const Geometry> a, const std::shared_ptr<const Matrix> b, const std::shared_ptr<const Matrix> c,
         const std::shared_ptr<const Matrix> ocoeff, const bool store = false, const bool rhf = false, const double scale_ex = 1.0)
#endif

  hcore_->print("hcore", 10);
  shared_ptr<Matrix> intermediate = make_shared<Matrix>(*tildex_ % *hcore_ * *tildex_);
  intermediate->diagonalize(eig_.get());
  shared_ptr<const Matrix> coeff = make_shared<Matrix>(*tildex_ * *intermediate);

  for (int i = 0; i != tildex_->mdim(); ++i)
    cout << eig_[i] << endl; 

  cout << nocc_ << endl;

  auto fock = make_shared<const Fock<1>>(geom_, hcore_, shared_ptr<const Matrix>(), coeff->slice(0, nocc_), false/*do_grad*/, true/*rhf*/);
  fock->print();



}
