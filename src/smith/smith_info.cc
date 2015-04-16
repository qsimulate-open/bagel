//
// BAGEL - Parallel electron correlation program.
// Filename: smith_info.cc
// Copyright (C) 2015 Toru Shiozaki
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

#include <src/smith/smith_info.h>
#include <src/multi/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;

template<>
tuple<shared_ptr<const RDM<1>>, shared_ptr<const RDM<2>>> SMITH_Info<double>::rdm12(const int ist, const int jst) const {
  return ref_->rdm12(ist, jst);
}


template<>
tuple<shared_ptr<const RDM<3>>, shared_ptr<const RDM<4>>> SMITH_Info<double>::rdm34(const int ist, const int jst) const {
  return ref_->rdm34(ist, jst);
}


template<>
tuple<shared_ptr<const ZRDM<1>>, shared_ptr<const ZRDM<2>>> SMITH_Info<complex<double>>::rdm12(const int ist, const int jst) const {
  assert(false);
  return make_tuple(shared_ptr<const ZRDM<1>>(), shared_ptr<const ZRDM<2>>());
}


template<>
tuple<shared_ptr<const ZRDM<3>>, shared_ptr<const ZRDM<4>>> SMITH_Info<complex<double>>::rdm34(const int ist, const int jst) const {
  assert(false);
  return make_tuple(shared_ptr<const ZRDM<3>>(), shared_ptr<const ZRDM<4>>());
}


template<>
shared_ptr<const RDM<1>> SMITH_Info<double>::rdm1_av() const {
  return ref_->rdm1_av();
}


template<>
shared_ptr<const ZRDM<1>> SMITH_Info<complex<double>>::rdm1_av() const {
  return nullptr;
}


template<>
shared_ptr<const Matrix> SMITH_Info<double>::coeff() const {
  return ref_->coeff();
}


template<>
shared_ptr<const ZMatrix> SMITH_Info<complex<double>>::coeff() const {
  shared_ptr<const ZMatrix> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
  // reformat so that it will be blocked
  return ZCASSCF::format_coeff(nclosed(), nact(), nvirt(), c, true);
}


template<>
shared_ptr<const Matrix> SMITH_Info<double>::hcore() const {
  return ref_->hcore();
}


template<>
shared_ptr<const ZMatrix> SMITH_Info<complex<double>>::hcore() const {
  // TODO implement
  assert(false);
  return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class SMITH_Info<double>;
template class SMITH_Info<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
