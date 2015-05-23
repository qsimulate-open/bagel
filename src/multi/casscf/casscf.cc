//
// BAGEL - Parallel electron correlation program.
// Filename: casscf.cc
// Copyright (C) 2011 Toru Shiozaki
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
#include <src/scf/hf/fock.h>
#include <src/multi/casscf/casscf.h>
#include <src/multi/casscf/qvec.h>
#include <src/wfn/construct_method.h>
#include <src/wfn/localization.h>

using namespace std;
using namespace bagel;


CASSCF::CASSCF(std::shared_ptr<const PTree> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : Method(idat, geom, re), hcore_(make_shared<Hcore>(geom)) {

  const shared_ptr<const PTree> mdata = idat->get_child_optional("molecule");
  if (ref_ && mdata) {
    //update geometry
    geom_ = make_shared<Geometry>(*geom, mdata);
    hcore_ = make_shared<Hcore>(geom_);
    //guess reference
    guess_ref_ = re->project_coeff(geom_, false);
    shared_ptr<const Matrix> projected = guess_ref_->coeff();
    // orthonormalize the "projected" coefficients
    Overlap S(geom_);
    Matrix S_invhalf = (*projected) % S * (*projected);
    S_invhalf.inverse_half();
    auto coeff = make_shared<Matrix>(*projected * S_invhalf);
    //update reference
    ref_ = make_shared<Reference>(geom_, make_shared<const Coeff>(move(*coeff)), guess_ref_->nclosed(), /*0*/guess_ref_->nact(), guess_ref_->nvirt());
    //do SCF
    auto hfdata = idat->get_child_optional("hf") ? idat->get_child_optional("hf") : make_shared<PTree>();
    auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", hfdata, geom_, ref_));
    rhf->compute();
    ref_ = rhf->conv_to_ref();
    //localization
    shared_ptr<const PTree> localize_data = idat->get_child_optional("localization");
    if (localize_data) {
      string localizemethod = localize_data->get<string>("algorithm", "pm");
      shared_ptr<OrbitalLocalization> localization;
      if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
        localization = make_shared<PMLocalization>(localize_data, ref_);
      else throw runtime_error("Unrecognized orbital localization method");
      //update reference with localized orbital
      shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize());
      ref_ = make_shared<const Reference>(*ref_, new_coeff);
    }
  }

  if (!ref_) {
    auto scf = make_shared<RHF>(idat, geom);
    scf->compute();
    ref_ = scf->conv_to_ref();
  }

  common_init();

}


void CASSCF::common_init() {
  // at the moment I only care about C1 symmetry, with dynamics in mind
  if (geom_->nirrep() > 1) throw runtime_error("CASSCF: C1 only at the moment.");
  print_header();

  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    set<int> active_indices;
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive) active_indices.insert(lexical_cast<int>(i->data()) - 1);

    if (guess_ref_) {
      const int nbasis = geom_->nbasis();
      //project to find active orbitals
      const int nclosed_hf = ref_->nclosed();
      const int nvirt_hf = ref_->nvirt();

      guess_ref_ = guess_ref_->set_active(active_indices);
      const int nclosed = guess_ref_->nclosed();
      const int nact = guess_ref_->nact();
      const int nvirt = nbasis - nclosed - nact;
      cout << nclosed_hf << " " << nvirt_hf << endl;
      cout << nclosed << " " << nact << " " << nvirt << endl;
      cout << nbasis << endl;
      assert(nbasis == nclosed_hf + nvirt_hf);
      assert(nbasis == nclosed + nact + nvirt);

      vector<tuple<shared_ptr<const Matrix>, pair<int, int>, int, string, bool>> svd_info;

      auto active = make_shared<Matrix>(nbasis, nact);
      active->copy_block(0, 0, nbasis, nact, guess_ref_->coeff()->get_submatrix(0, nclosed, nbasis, nact));
      svd_info.emplace_back(active, make_pair(0, nclosed_hf), nclosed_hf - nclosed, "closed subspace", true);
      svd_info.emplace_back(active, make_pair(nclosed_hf, nbasis), nvirt_hf - nvirt, "virtual subspace", false);
      assert(nclosed_hf - nclosed + nvirt_hf - nvirt == nact);

      Overlap S(geom_);

      shared_ptr<Matrix> out_coeff = ref_->coeff()->copy();
      size_t closed_position = 0;
      size_t active_position = nclosed;
      size_t virt_position = nclosed + nact;

      for (auto& subset : svd_info) {
        const Matrix& active = *get<0>(subset);
        pair<int, int> bounds = get<1>(subset);
        const int norb = get<2>(subset);
        const string set_name = get<3>(subset);
        const bool closed = get<4>(subset);

        shared_ptr<Matrix> subcoeff = ref_->coeff()->slice_copy(bounds.first, bounds.second);

        const Matrix overlaps(active % S * *subcoeff);

        multimap<double, int> norms;

        for(int i = 0; i < overlaps.mdim(); ++i) {
          const double norm = blas::dot_product(overlaps.element_ptr(0, i), overlaps.ndim(), overlaps.element_ptr(0, i));
          norms.emplace(norm, i);
        }

        active_thresh_ = idata_->get<double>("active_thresh", 0.5);
        cout << endl << "  o Forming CASSCF active orbitals arising from " << (closed ? "closed " : "virtual ") << "guess orbitals. Threshold for inclusion in cadidate space: " << setw(6) << setprecision(3) << active_thresh_ << endl;

        vector<int> active_list;
        double max_overlap, min_overlap;
        {
          auto end = norms.rbegin(); advance(end, norb);
          end = find_if(end, norms.rend(), [this] (const pair<const double, int>& p) { return p.first < active_thresh_; });
          for_each(norms.rbegin(), end, [&active_list] (const pair<const double, int>& p) { active_list.emplace_back(p.second); });
          auto mnmx = minmax_element(norms.rbegin(), end);
          tie(min_overlap, max_overlap) = make_tuple(mnmx.first->first, mnmx.second->first);
        }

        const int active_size = active_list.size();
        cout << "    - size of candidate space: " << active_size << endl;
        cout << "    - largest overlap with guess space: " << max_overlap << ", smallest: " << min_overlap << endl;

        if (active_size != norb) {
          throw runtime_error("Try adjust active_thresh.");
        }
        else {
          set<int> active_set(active_list.begin(), active_list.end());
          for (size_t i = 0; i < subcoeff->mdim(); ++i)
            if (active_set.count(i) == 0)
              copy_n(subcoeff->element_ptr(0, i), nbasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
            else
              copy_n(subcoeff->element_ptr(0, i), nbasis, out_coeff->element_ptr(0, active_position++));
        }
      }
      //Reordered reference
      ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nvirt);
    }
    else {
      ref_ = ref_->set_active(active_indices);
      cout << " " << endl;
      cout << "    ==== Active orbitals : ===== " << endl;
      for (auto& i : active_indices) cout << "         Orbital " << i+1 << endl;
      cout << "    ============================ " << endl << endl;
    }
  }

  // first set coefficient
  coeff_ = ref_->coeff();
  if (geom_->nbasis() != coeff_->mdim()) {
    Overlap ovl(geom_);
    shared_ptr<const Matrix> tildex = ovl.tildex();

    Matrix c(coeff_->ndim(), tildex->mdim());
    c.copy_block(0, 0, coeff_->ndim(), coeff_->mdim(), coeff_);

    shared_ptr<const Matrix> trans = get<0>((*tildex % ovl * *coeff_).svd());
    c.copy_block(0, coeff_->mdim(), coeff_->ndim(), tildex->mdim()-coeff_->mdim(), *tildex * trans->slice(coeff_->mdim(), tildex->mdim()));
    coeff_ = make_shared<Coeff>(move(c));

#ifndef NDEBUG
    Matrix unit(coeff_->mdim(), coeff_->mdim()); unit.unit();
    assert((*coeff_ % ovl * *coeff_ - unit).rms() < 1.0e-10);
#endif
  }

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 50);
  // get maxiter from the input
  max_micro_iter_ = idata_->get<int>("maxiter_micro", 100);
  // get nstate from the input
  nstate_ = idata_->get<int>("nstate", 1);
  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", 5.0e-6);
  // option for printing natural orbitals
  natocc_ = idata_->get<bool>("natocc",false);

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

  nbasis_ = coeff_->mdim();
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  const int idel = geom_->nbasis() - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;


  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and total orbitals.
  mute_stdcout();
  if (nact_) {
    auto idata = make_shared<PTree>(*idata_);
    idata->erase("active");
    fci_ = make_shared<KnowlesHandy>(idata, geom_, ref_, nclosed_, nact_); // nstate does not need to be specified as it is in idata_...
  }
  resume_stdcout();


  schwarz_ = geom_->schwarz();

  cout <<  "  === CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


CASSCF::~CASSCF() {

}


void CASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}

void CASSCF::print_iteration(int iter, int miter, int tcount, const vector<double> energy, const double error, const double time) const {
  if (energy.size() != 1 && iter) cout << endl;

  int i = 0;
  for (auto& e : energy) {
    cout << "  " << setw(5) << iter << setw(3) << i << setw(4) << miter << setw(4) << tcount
                 << setw(16) << fixed << setprecision(8) << e << "   "
                 << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2)
                 << time << endl;
    ++i;
  }
}

static streambuf* backup_stream_;
static ofstream* ofs_;

void CASSCF::mute_stdcout() {
  ofstream* ofs(new ofstream("casscf.log",(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
}


void CASSCF::resume_stdcout() {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


shared_ptr<Matrix> CASSCF::ao_rdm1(shared_ptr<const RDM<1>> rdm1, const bool inactive_only) const {
  // first make 1RDM in MO
  const size_t nmobasis = coeff_->mdim();
  auto mo_rdm1 = make_shared<Matrix>(nmobasis, nmobasis);
  for (int i = 0; i != nclosed_; ++i) mo_rdm1->element(i,i) = 2.0;
  if (!inactive_only) {
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nact_; ++j) {
        mo_rdm1->element(nclosed_+j, nclosed_+i) = rdm1->element(j,i);
      }
    }
  }
  // transform into AO basis
  return make_shared<Matrix>(*coeff_ * *mo_rdm1 ^ *coeff_);
}



void CASSCF::one_body_operators(shared_ptr<Matrix>& f, shared_ptr<Matrix>& fact, shared_ptr<Matrix>& factp, shared_ptr<Matrix>& gaa,
                                shared_ptr<RotFile>& d, const bool superci) const {

  shared_ptr<Matrix> finact;

  // get quantity Q_xr = 2(xs|tu)P_rs,tu (x=general)
  // note: this should be after natorb transformation.
  auto qxr = make_shared<Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

  {
    // Fock operators
    // make a matrix that contains rdm1_av
    auto rdm1mat = make_shared<Matrix>(nact_, nact_);
    copy_n(fci_->rdm1_av()->data(), rdm1mat->size(), rdm1mat->data());
    rdm1mat->sqrt();
    rdm1mat->scale(1.0/sqrt(2.0));
    auto acoeff = coeff_->slice(nclosed_, nclosed_+nact_);

    finact = make_shared<Matrix>(*coeff_ % *fci_->jop()->core_fock() * *coeff_);
    auto fact_ao = make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeff * *rdm1mat, false, /*rhf*/true);
    f = make_shared<Matrix>(*finact + *coeff_% *fact_ao * *coeff_);
  }
  {
    // active-x Fock operator Dts finact_sx + Qtx
    fact = qxr->copy();// nbasis_ runs first
    for (int i = 0; i != nact_; ++i)
      daxpy_(nbasis_, occup_(i), finact->element_ptr(0,nclosed_+i), 1, fact->data()+i*nbasis_, 1);
  }

  {
    // active Fock' operator (Fts+Fst) / (ns+nt)
    factp = make_shared<Matrix>(nact_, nact_);
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nact_; ++j) {
#if 1
        if (occup_(i)+occup_(j) > occup_thresh)
          factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) / (occup_(i)+occup_(j));
        else
          factp->element(j,i) = 0.0;
#else
        factp->element(j,i) = (fact->element(j+nclosed_,i)+fact->element(i+nclosed_,j)) *0.5;
#endif
      }
  }

  // G matrix (active-active) Drs,tu Factp_tu - delta_rs nr sum_v Factp_vv
  gaa = factp->clone();
  dgemv_("N", nact_*nact_, nact_*nact_, 1.0, fci_->rdm2_av()->data(), nact_*nact_, factp->data(), 1, 0.0, gaa->data(), 1);
  double p = 0.0;
  for (int i = 0; i != nact_; ++i) p += occup_(i) * factp->element(i,i);
  for (int i = 0; i != nact_; ++i) gaa->element(i,i) -= occup_(i) * p;

  // denominator
  auto denom = make_shared<RotFile>(nclosed_, nact_, nvirt_);
  fill_n(denom->data(), denom->size(), 1.0e100);

  double* target = denom->ptr_va();
  for (int i = 0; i != nact_; ++i) {
    if (occup_(i) > occup_thresh) {
      for (int j = 0; j != nvirt_; ++j, ++target)
        *target = (gaa->element(i,i) + occup_(i)*f->element(j+nocc_, j+nocc_)) / (superci ? occup_(i) : 1.0);
    } else {
      for (int j = 0; j != nvirt_; ++j, ++target)
        *target = 1.0/occup_thresh;
    }
  }

  target = denom->ptr_vc();
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nvirt_; ++j, ++target)
      *target = (f->element(j+nocc_, j+nocc_) - f->element(i, i)) / (superci ? 2.0 : 1.0);

  target = denom->ptr_ca();
  for (int i = 0; i != nact_; ++i) {
    if (2.0-occup_(i) > occup_thresh) {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = ((f->element(nclosed_+i,nclosed_+i)*2.0-fact->element(i+nclosed_,i)) - f->element(j, j)*(2.0-occup_(i))) / (superci ? 2.0-occup_(i) : 1.0);
    } else {
      for (int j = 0; j != nclosed_; ++j, ++target)
        *target = 1.0/occup_thresh;
    }
  }
  d = denom;
}


shared_ptr<const Coeff> CASSCF::update_coeff(const shared_ptr<const Matrix> cold, shared_ptr<const Matrix> mat) const {
  auto cnew = make_shared<Coeff>(*cold);
  int nbas = cold->ndim();
  assert(nbas == geom_->nbasis());
  dgemm_("N", "N", nbas, nact_, nact_, 1.0, cold->data()+nbas*nclosed_, nbas, mat->data(), nact_,
                   0.0, cnew->data()+nbas*nclosed_, nbas);
  return cnew;
}


shared_ptr<const Coeff> CASSCF::semi_canonical_orb() const {
  auto rdm1mat = make_shared<Matrix>(nact_, nact_);
  copy_n(fci_->rdm1_av()->data(), rdm1mat->size(), rdm1mat->data());
  rdm1mat->sqrt();
  rdm1mat->scale(1.0/sqrt(2.0));
  auto ocoeff = coeff_->slice(0, nclosed_);
  auto acoeff = coeff_->slice(nclosed_, nocc_);
  auto vcoeff = coeff_->slice(nocc_, nbasis_);

  VectorB eig(coeff_->mdim());
  Fock<1> fock(geom_, fci_->jop()->core_fock(), nullptr, acoeff * *rdm1mat, false, /*rhf*/true);
  Matrix trans(nbasis_, nbasis_);
  trans.unit();
  if (nclosed_) {
    Matrix ofock = ocoeff % fock * ocoeff;
    ofock.diagonalize(eig);
    trans.copy_block(0, 0, nclosed_, nclosed_, ofock);
  }
  Matrix vfock = vcoeff % fock * vcoeff;
  vfock.diagonalize(eig);
  trans.copy_block(nocc_, nocc_, nvirt_, nvirt_, vfock);
  return make_shared<Coeff>(*coeff_ * trans);
}


shared_ptr<const Reference> CASSCF::conv_to_ref() const {
  shared_ptr<Reference> out;
  if (nact_) {
    out = make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_av(),
                                 fci_->rdm1(), fci_->rdm2(), fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
    // TODO
    // compute one-body operators
    shared_ptr<Matrix> f;
    shared_ptr<Matrix> fact, factp, gaa;
    shared_ptr<RotFile>  denom;
    one_body_operators(f, fact, factp, gaa, denom);
    if (natocc_) print_natocc();

    *f *= 2.0;

    for (int i = 0; i != nbasis_; ++i) {
      for (int j = 0; j != nbasis_; ++j) {
        if (i < nocc_ && j < nocc_) continue;
        f->element(j,i) = 0.0;
      }
    }
    for (int j = 0; j != nact_; ++j) {
      for (int i = 0; i != nocc_; ++i) {
        f->element(i,j+nclosed_) = fact->element(i,j);
      }
    }

    auto erdm = make_shared<Matrix>(*coeff_ * *f ^ *coeff_);

    out->set_erdm1(erdm);
    out->set_nstate(nstate_);
  } else {
    out = make_shared<Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_av());
  }
  return out;
}


void CASSCF::print_natocc() const {
  assert(occup_.size() > 0);
  cout << " " << endl;
  cout << "  ========       state-averaged       ======== " << endl;
  cout << "  ======== natural occupation numbers ======== " << endl;
  for (int i=0; i!=occup_.size(); ++i)
    cout << setprecision(4) << "   Orbital " << i << " : " << occup_[i] << endl;
  cout << "  ============================================ " << endl;
}
