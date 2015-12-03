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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/smith_info.h>
#include <src/wfn/relcoeff.h>

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
tuple<shared_ptr<const RDM<3>>, shared_ptr<const RDM<3>>> SMITH_Info<double>::rdm34f(const int ist, const int jst, shared_ptr<const Matrix> fock) const {
  return ref_->rdm34f(ist, jst, fock);
}


template<>
tuple<shared_ptr<const ZRDM<1>>, shared_ptr<const ZRDM<2>>>
  SMITH_Info<complex<double>>::rdm12(const int ist, const int jst) const {

  auto ref = dynamic_pointer_cast<const RelReference>(ref_);
  auto rdm1 = ref->rdm1(ist, jst);
  auto rdm2 = ref->rdm2(ist, jst);
  return make_tuple(expand_kramers(rdm1,nact()), expand_kramers(rdm2,nact()));
}


template<>
tuple<shared_ptr<const ZRDM<3>>, shared_ptr<const ZRDM<4>>>
  SMITH_Info<complex<double>>::rdm34(const int ist, const int jst) const {

  auto ref = dynamic_pointer_cast<const RelReference>(ref_);
  auto rdm3 = ref->rdm3(ist, jst);
  auto rdm4 = ref->rdm4(ist, jst);
  return make_tuple(expand_kramers(rdm3,nact()), expand_kramers(rdm4,nact()));
}


template<>
tuple<shared_ptr<const ZRDM<3>>, shared_ptr<const ZRDM<3>>>
  SMITH_Info<complex<double>>::rdm34f(const int ist, const int jst, shared_ptr<const ZMatrix> fock) const {

  // TODO not the best code
  shared_ptr<const ZRDM<3>> rdm3;
  shared_ptr<const ZRDM<4>> rdm4;
  tie(rdm3, rdm4) = rdm34(ist, jst);

  shared_ptr<ZRDM<3>> rdm4f = rdm3->clone();;
  auto rdm4v = group(group(*rdm4, 6,8), 0,6);
  auto rdm4vf = group(*rdm4f, 0, 6);
  contract(1.0, rdm4v, {0,1}, group(*fock,0,2), {1}, 0.0, rdm4vf, {0});

  return make_tuple(rdm3, rdm4f);
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
shared_ptr<const CIWfn> SMITH_Info<double>::ciwfn() const {
  return ref_->ciwfn();
}


template<>
shared_ptr<const RelCIWfn> SMITH_Info<complex<double>>::ciwfn() const {
  return dynamic_pointer_cast<const RelReference>(ref_)->ciwfn();
}


template<>
shared_ptr<const Matrix> SMITH_Info<double>::coeff() const {
  return ref_->coeff();
}


template<>
shared_ptr<const ZMatrix> SMITH_Info<complex<double>>::coeff() const {
  shared_ptr<const RelCoeff_Striped> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
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

#endif
