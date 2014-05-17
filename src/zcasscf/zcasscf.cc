//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <fstream>
#include <src/rel/dirac.h>
#include <src/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;

ZCASSCF::ZCASSCF(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref)
  : Method(idat, geom, ref) {

  if (ref != nullptr && ref->coeff()->ndim() == geom->nbasis()) {
    nr_coeff_ = ref->coeff();
  }
// relref needed for many things below ; TODO eliminate dependence on ref_ being a relref
  {
    auto idata_tmp = make_shared<PTree>(*idata_);
    const int ctmp = idata_->get<int>("charge", 0);
    const int nele = geom->nele();
    if ((nele - ctmp) % 2) {
      idata_tmp->erase("charge");
      idata_tmp->put<int>("charge", ctmp - 1);
    }
    auto scf = make_shared<Dirac>(idata_tmp, geom_);
    scf->compute();
    ref_ = scf->conv_to_ref();
  }

  init();

}


void ZCASSCF::init() {
  print_header();

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  if (!geom_->dfs())
    geom_ = geom_->relativistic(relref->gaunt());

  nneg_ = relref->nneg();

  // first set coefficient
  {
    shared_ptr<const ZMatrix> ctmp = relref->relcoeff_full();
    shared_ptr<ZMatrix> coeff = ctmp->clone();
    const int npos = ctmp->mdim() - nneg_;
    coeff->copy_block(0, 0, ctmp->mdim(), npos, ctmp->slice(nneg_, nneg_+npos));
    coeff->copy_block(0, npos, ctmp->mdim(), nneg_, ctmp->slice(0, nneg_));
    coeff_ = coeff;
  }

  // set hcore and overlap
  hcore_   = make_shared<RelHcore>(geom_);
  overlap_ = make_shared<RelOverlap>(geom_);

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 100);
  // get maxiter from the input
  max_micro_iter_ = idata_->get<int>("maxiter_micro", 100);
  // get nstate from the input
  nstate_ = idata_->get<int>("nstate", 1);
#if 0
  // get istate from the input (for geometry optimization)
  istate_ = idata_->get<int>("istate", 0);
#endif
  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", thresh_);

  // nocc from the input. If not present, full valence active space is generated.
  nact_ = idata_->get<int>("nact", 0);
  nact_ = idata_->get<int>("nact_cas", nact_);

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = idata_->get<int>("nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore_only() / 2;
  }
  nocc_ = nclosed_ + nact_;

  nbasis_ = coeff_->mdim()/2;
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");
  nvirtnr_ = nvirt_ - nneg_/2;

  charge_ = idata_->get<int>("charge", 0);
  if (nclosed_*2 > geom_->nele() - charge_)
    throw runtime_error("two many closed orbitals in the input");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  gaunt_ = relref->gaunt();
  breit_ = relref->breit();
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;

  const int idel = geom_->nbasis()*2 - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;

  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and active orbitals.
  if (nact_) {
    mute_stdcout(/*fci*/true);
    fci_ = make_shared<ZHarrison>(idata_, geom_, ref_, nclosed_, nact_); // nstate does not need to be specified as it is in idata_...
    resume_stdcout();
  }

  cout <<  "  === Dirac CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


void ZCASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}


void ZCASSCF::print_iteration(int iter, int miter, int tcount, const vector<double> energy, const double error, const double time) const {
  if (energy.size() != 1 && iter) cout << endl;
  int i = 0;
  cout << "Cycle" << setw(5) << iter << setw(3) << i << setw(4) << miter << setw(4) << tcount
               << setw(20) << fixed << setprecision(12) << energy[(energy.size() > 0 ? energy.size()-1 : 0)] << "   "
               << setw(10) << scientific << setprecision(4) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2)
               << time << endl;
}


static streambuf* backup_stream_;
static ofstream* ofs_;

void ZCASSCF::mute_stdcout(const bool fci) const {
  if (fci) {
    ofstream* ofs(new ofstream("casscf.log",(backup_stream_ ? ios::app : ios::trunc)));
    ofs_ = ofs;
    backup_stream_ = cout.rdbuf(ofs->rdbuf());
  } else {
    ofstream* ofs(new ofstream("microiter.log",(backup_stream_ ? ios::app : ios::trunc)));
    ofs_ = ofs;
    backup_stream_ = cout.rdbuf(ofs->rdbuf());
  }
}


void ZCASSCF::resume_stdcout() const {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


shared_ptr<const ZMatrix> ZCASSCF::transform_rdm1() const {
  assert(fci_);
  shared_ptr<const ZMatrix> out = fci_->rdm1_av();
  assert(out->ndim() == 2*nact_ && out->mdim() == 2*nact_);
  return out;
}


shared_ptr<const ZMatrix> ZCASSCF::active_fock(shared_ptr<const ZMatrix> rdm1) const {
  // form natural orbitals
  unique_ptr<double[]> eig(new double[nact_*2]);
  auto tmp = make_shared<ZMatrix>(*rdm1);
  tmp->diagonalize(eig.get());
  auto ocoeff = coeff_->slice(nclosed_*2, nclosed_*2+nact_*2);
  // D_rs = C*_ri D_ij (C*_rj)^+. Dij = U_ik L_k (U_jk)^+. So, C'_ri = C_ri * U*_ik
  auto natorb = make_shared<ZMatrix>(*ocoeff * *tmp->get_conjg());

  // scale using eigen values
  for (int i = 0; i != nact_*2; ++i) {
    assert(eig[i] >= -1.0e-14);
    const double fac = eig[i] > 0 ? sqrt(eig[i]) : 0.0;
    for_each(natorb->element_ptr(0, i), natorb->element_ptr(0, i+1), [&fac](complex<double>& a) { a *= fac; });
  }

  auto zero = make_shared<ZMatrix>(geom_->nbasis()*4, geom_->nbasis()*4);
  return make_shared<const DFock>(geom_, zero, natorb, gaunt_, breit_, /*store half*/false, /*robust*/breit_);
}


shared_ptr<ZMatrix> ZCASSCF::make_natural_orbitals(shared_ptr<const ZMatrix> rdm1) {
  // input should be 1rdm in kramers format
  shared_ptr<ZMatrix> tmp = rdm1->copy();

  unique_ptr<double[]> vec(new double[rdm1->ndim()]);
  zquatev_(tmp->ndim(), tmp->data(), vec.get()); // TODO : maybe replace with standard diagonalize

  map<int,int> emap;
  auto buf2 = tmp->clone();
  vector<double> vec2(tmp->ndim());
  // sort eigenvectors so that buf is close to a unit matrix
  // target column
  for (int i = 0; i != tmp->ndim(); ++i) {
    // first find the source column
    tuple<int, double> max = make_tuple(-1, 0.0);
    for (int j = 0; j != tmp->ndim(); ++j)
      if (sqrt(real(tmp->element(i,j)*tmp->get_conjg()->element(i,j))) > get<1>(max))
        max = make_tuple(j, sqrt(real(tmp->element(i,j)*tmp->get_conjg()->element(i,j))));

    // register to emap
    if (emap.find(get<0>(max)) != emap.end()) throw logic_error("this should not happen. make_natural_orbitals()");
    emap.insert(make_pair(get<0>(max), i));

    // copy to the target
    copy_n(tmp->element_ptr(0,get<0>(max)), tmp->ndim(), buf2->element_ptr(0,i));
    vec2[i] = vec[get<0>(max)];
  }

  // fix the phase
  for (int i = 0; i != tmp->ndim(); ++i) {
    if (real(buf2->element(i,i)) < 0.0)
      blas::scale_n(-1.0, buf2->element_ptr(0,i), tmp->ndim());
  }
  // copy eigenvalues TODO: change to blas
  for (int i=0; i!=tmp->ndim()/2; ++i)
    vec2[tmp->ndim()/2 + i] = vec2[i];
  occup_ = vec2;
  coeff_ = update_coeff(coeff_, buf2);

  return buf2;
}


shared_ptr<const ZMatrix> ZCASSCF::natorb_rdm1_transform(const shared_ptr<ZMatrix> coeff, shared_ptr<const ZMatrix> rdm1) const {
  shared_ptr<ZMatrix> tmp = rdm1->clone();
  const complex<double>* start = coeff->data();
  int ndim = coeff->ndim();
  unique_ptr<complex<double>[]> buf(new complex<double>[ndim*ndim]);
  zgemm3m_("N", "N", ndim, ndim, ndim, 1.0, rdm1->data(), ndim, start, ndim, 0.0, buf.get(), ndim);
  zgemm3m_("C", "N", ndim, ndim, ndim, 1.0, start, ndim, buf.get(), ndim, 0.0, tmp->data(), ndim);
  auto out = make_shared<const ZMatrix>(*tmp);
  return out;
}


shared_ptr<const ZMatrix> ZCASSCF::natorb_rdm2_transform(const shared_ptr<ZMatrix> coeff, shared_ptr<const ZMatrix> rdm2) const {
  shared_ptr<ZMatrix> tmp = rdm2->clone();
  const complex<double>* start = coeff->data();
  int ndim  = coeff->ndim();
  int ndim2 = rdm2->ndim();
  unique_ptr<complex<double>[]> buf(new complex<double>[ndim2*ndim2]);
  // first half transformation 
  zgemm3m_("N", "N", ndim2*ndim, ndim, ndim, 1.0, rdm2->data(), ndim2*ndim, start, ndim, 0.0, buf.get(), ndim2*ndim);
  for (int i = 0; i != ndim; ++i)
    zgemm3m_("N", "C", ndim2, ndim, ndim, 1.0, buf.get()+i*ndim2*ndim, ndim2, start, ndim, 0.0, tmp->data()+i*ndim2*ndim, ndim2);
  // then tranpose
  blas::transpose(tmp->data(), ndim2, ndim2, buf.get());
  // and do it again
  zgemm3m_("N", "N", ndim2*ndim, ndim, ndim, 1.0, buf.get(), ndim2*ndim, start, ndim, 0.0, tmp->data(), ndim2*ndim);
  for (int i = 0; i != ndim; ++i)
    zgemm3m_("N", "C", ndim2, ndim, ndim, 1.0, tmp->data()+i*ndim2*ndim, ndim2, start, ndim, 0.0, buf.get()+i*ndim2*ndim, ndim2);
  // to make sure for non-symmetric density matrices (and anyway this should be cheap).
  blas::transpose(buf.get(), ndim2, ndim2, tmp->data());
  auto out = make_shared<const ZMatrix>(*tmp);
  return out;
}


shared_ptr<const ZMatrix> ZCASSCF::update_coeff(shared_ptr<const ZMatrix> cold, shared_ptr<const ZMatrix> natorb) const {
  // see active_fock for explanation of conjugation for natorb
  auto cnew = make_shared<ZMatrix>(*cold);
  int n    = natorb->ndim();
  int nbas = cold->ndim();
  zgemm3m_("N", "N", nbas, n, n, 1.0, cold->data()+nbas*nclosed_*2, nbas, natorb->get_conjg()->data(), n,
                   0.0, cnew->data()+nbas*nclosed_*2, nbas);
  return cnew;
}


shared_ptr<const ZMatrix> ZCASSCF::update_qvec(shared_ptr<const ZMatrix> qold, shared_ptr<const ZMatrix> natorb) const {
  auto qnew = make_shared<ZMatrix>(*qold);
  int n    = natorb->ndim();
  int nbas = qold->ndim();
  // first transformation
  zgemm3m_("N", "N", nbas, n, n, 1.0, qold->data(), nbas, natorb->data(), n, 0.0, qnew->data(), nbas);
  // second transformation for the active-active block
  auto qtmp = qnew->get_submatrix(nclosed_*2, 0, n, n)->copy(); 
  *qtmp = *natorb % *qtmp;
  qnew->copy_block(nclosed_*2, 0, n, n, qtmp->data());
  return qnew;
}
