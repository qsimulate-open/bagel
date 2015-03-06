//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/casscf.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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
#include <src/scf/hf/fock.h>
#include <src/asd/orbital/oo.h>

#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

ASD_OO::ASD_OO(shared_ptr<const PTree> idat, shared_ptr<Dimer> dimer)
  : Method(idat, dimer->sgeom(), dimer->sref()), dimer_(dimer), hcore_(make_shared<Hcore>(dimer->sgeom())) {
  common_init();
}


void ASD_OO::common_init() {
  print_header();

  // first set coefficient
  coeff_ = ref_->coeff();

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 50);
  // get nstate from the ASD input
  nstate_ = idata_->get_child_optional("asd")->get<int>("nstates", 1);
  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh
  precond_ = idata_->get<double>("precondition", 5.0e-4);

  // active space
  nact_ = ref_->nact();
  nactA_ = dimer_->active_refs().first->nact();
  nactB_ = dimer_->active_refs().second->nact();
  rasA_ = {0, nactA_, 0};
  rasB_ = {0, nactB_, 0};
  if (idata_->get_child_optional("asd")->get<string>("method") == "ras") {
    //TODO
  }

  assert(nactA_ + nactB_ == nact_);

  nclosed_ = ref_->nclosed();

  nocc_ = nclosed_ + nact_;

  nbasis_ = coeff_->mdim();
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    *  unit A  : " << setw(6) << nactA_ << endl;
  cout << "    *  unit B  : " << setw(6) << nactB_ << endl;
  cout << "    * nocc     : " << setw(6) << nocc_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;

  cout <<  "  === ASD Orbital Optimization iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


ASD_OO::~ASD_OO() {

}


void ASD_OO::print_header() const {
  cout << "  --------------------------------------------" << endl;
  cout << "      ASD Orbital Optimization calculation    " << endl;
  cout << "  --------------------------------------------" << endl << endl;
}

void ASD_OO::print_iteration(int iter, int miter, int tcount, const vector<double> energy, const double error, const double time) const {
  if (energy.size() != 1 && iter) cout << endl;

  int i = 0;
  for (auto& e : energy) {
    cout << "  " << setw(5) << iter << setw(3) << i << setw(4) << miter << setw(4) << tcount
                 << setw(20) << fixed << setprecision(12) << e << "   "
                 << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2)
                 << time << endl;
    ++i;
  }
}


static streambuf* backup_stream_;
static ofstream* ofs_;


void ASD_OO::mute_stdcout() {
  ofstream* ofs(new ofstream("asd_casscf.log",(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
}


void ASD_OO::resume_stdcout() {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


shared_ptr<Matrix> ASD_OO::Qvec(const int n, const int m, shared_ptr<const Matrix> coeff, const size_t nclosed) const {
  assert(n == coeff->mdim());

  std::shared_ptr<DFHalfDist> half;
  const MatView cdata = coeff_->slice(nclosed_, nclosed_+nact_);
  half = geom_->df()->compute_half_transform(cdata);

  // J^{-1}(D|xy)
  // TODO : DFDistT needs to be modified to handle cases where number of nodes is larger than half->nocc() * cdata.mdim()
  shared_ptr<const DFFullDist> full;
  if (half->nocc() * coeff->mdim() > mpi__->size()) {
    full = half->apply_JJ()->compute_second_transform(coeff->slice(nclosed, nclosed+m));
  } else {
    full = half->compute_second_transform(coeff->slice(nclosed, nclosed+m))->apply_JJ();
  }

  // [D|tu] = (D|xy)Gamma_xy,tu
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(*rdm2_);

  // (r,u) = (rt|D)[D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  // MO transformation of the first index
  auto out = make_shared<Matrix>(*coeff % *tmp);
  assert(n == out->ndim() && m == out->mdim());
  return out;

}


double ASD_OO::check_symmetric(std::shared_ptr<Matrix>& mat) const {
  int n = mat->ndim();
  assert(n == mat->mdim());
  auto tran = make_shared<Matrix>(n,n);
  tran = mat->transpose();
  auto subt = make_shared<Matrix>(n,n);
  *subt = *mat - *tran;
  return subt->rms();
}


shared_ptr<const Coeff> ASD_OO::update_coeff(const shared_ptr<const Matrix> cold, shared_ptr<const Matrix> mat) const {
  auto cnew = make_shared<Coeff>(*cold);
  int nbas = cold->ndim();
  assert(nbas == geom_->nbasis());
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, mat->data(), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  return cnew;
}


shared_ptr<const Reference> ASD_OO::conv_to_ref() const {
  cout << "ASD_OO: conv_to_ref: not yet implemented" << endl;
  //TODO
  assert(false);
}
