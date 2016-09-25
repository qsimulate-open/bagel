//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasscf.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <fstream>
#include <src/scf/dhf/dirac.h>
#include <src/scf/dhf/dfock.h>
#include <src/util/math/quatmatrix.h>
#include <src/multi/zcasscf/zcasscf.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/reloverlap_london.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

ZCASSCF::ZCASSCF(const shared_ptr<const PTree> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> ref)
  : Method(idat, geom, ref) {
  if (!dynamic_pointer_cast<const RelReference>(ref)) {
    if (ref != nullptr && ref->coeff()->ndim() == geom->nbasis()) {
      nr_coeff_ = ref->coeff();
      if (geom->magnetism()) // TODO Implement nr_coeff_ for GIAO basis sets
        throw runtime_error("So far only relativistic input coefficients can be used for ZCASSCF with magnetic field");
    }
  }

  // relref needed for many things below ; TODO eliminate dependence on ref_ being a relref
  if (!dynamic_pointer_cast<const RelReference>(ref)) {
    auto idata_tmp = make_shared<PTree>(*idata_);
    const int ctmp = idata_->get<int>("charge", 0);
    const int nele = geom->nele();
    if ((nele - ctmp) % 2) {
      idata_tmp->erase("charge");
      idata_tmp->put<int>("charge", ctmp - 1);
    }

    // With a CASSCF initial guess, the DHF wavefunction will not be used...
    if (nr_coeff_)
      idata_tmp->put<double>("thresh", 1.0e-4);

    auto scf = make_shared<Dirac>(idata_tmp, geom_, ref);
    scf->compute();
    ref_ = scf->conv_to_ref();
  }

  init();

}


void ZCASSCF::init() {
  print_header();

  auto relref = dynamic_pointer_cast<const RelReference>(ref_);

  gaunt_ = idata_->get<bool>("gaunt",relref->gaunt());
  breit_ = idata_->get<bool>("breit",relref->breit());

  if (!geom_->dfs() || (gaunt_ != relref->gaunt()))
    geom_ = geom_->relativistic(gaunt_);

  // Invoke Kramer's symmetry for any case without magnetic field
  tsymm_ = !geom_->magnetism();

  // coefficient parameters
  const bool kramers_coeff = idata_->get<bool>("kramers_coeff", relref->kramers());

  nneg_ = geom_->nbasis()*2;

  // set hcore and overlap
  if (!geom_->magnetism()) {
    hcore_   = make_shared<RelHcore>(geom_);
    overlap_ = make_shared<RelOverlap>(geom_);
  } else {
    hcore_ = make_shared<RelHcore_London>(geom_);
    overlap_ = make_shared<RelOverlap_London>(geom_);
  }

  // nocc from the input. If not present, full valence active space is generated.
  nact_ = idata_->get<int>("nact", 0);
  nact_ = idata_->get<int>("nact_cas", nact_);
  if (!nact_) energy_.resize(1);
  // option for printing natural orbital occupation numbers
  natocc_ = idata_->get<bool>("natocc",false);

  // nclosed from the input. If not present, full core space is generated.
  nclosed_ = idata_->get<int>("nclosed", -1);
  if (nclosed_ < -1) {
    throw runtime_error("It appears that nclosed < 0. Check nocc value.");
  } else if (nclosed_ == -1) {
    cout << "    * full core space generated for nclosed." << endl;
    nclosed_ = geom_->num_count_ncore_only() / 2;
  }
  nocc_ = nclosed_ + nact_;

  nbasis_ = geom_->nbasis()*2;
  nvirt_ = nbasis_ - nocc_;
  if (nvirt_ < 0) throw runtime_error("It appears that nvirt < 0. Check the nocc value");
  nvirtnr_ = nvirt_ - nneg_/2;

  charge_ = idata_->get<int>("charge", 0);
  if (nclosed_*2 > geom_->nele() - charge_)
    throw runtime_error("two many closed orbitals in the input");

  // set coefficient
  const bool hcore_guess = idata_->get<bool>("hcore_guess", false);
  shared_ptr<const RelCoeff_Striped> scoeff;
  if (hcore_guess) {
    auto s12 = overlap_->tildex(1.0e-10);
    auto hctmp = make_shared<ZMatrix>(*s12 % *hcore_ * *s12);
    VectorB eig(hctmp->ndim());
    hctmp->diagonalize(eig);
    scoeff = make_shared<const RelCoeff_Striped>(*s12 * *hctmp, nclosed_, nact_, nvirtnr_, nneg_, /*move_neg*/true);
  } else if (nr_coeff_ == nullptr) {
    scoeff = make_shared<const RelCoeff_Striped>(*relref->relcoeff_full(), nclosed_, nact_, nvirtnr_, nneg_);
  } else {
    scoeff = nonrel_to_relcoeff(nr_coeff_)->striped_format();
  }

  // get maxiter from the input
  max_iter_ = idata_->get<int>("maxiter", 100);
  // get maxiter from the input
  max_micro_iter_ = idata_->get<int>("maxiter_micro", 20);

  // get thresh (for macro iteration) from the input
  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // get thresh (for micro iteration) from the input
  thresh_micro_ = idata_->get<double>("thresh_micro", 5.0e-6);

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl;

  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;
  cout << "    * active space: " << geom_->nele() - charge_ - nclosed_*2 << " electrons in " << nact_ << " orbitals" << endl;
  cout << "    * time-reversal symmetry " << (tsymm_ ? "will be assumed." : "violation will be permitted.") << endl;

  const int idel = geom_->nbasis()*2 - nbasis_;
  if (idel)
    cout << "      Due to linear dependency, " << idel << (idel==1 ? " function is" : " functions are") << " omitted" << endl;

  // initialize coefficient to enforce kramers symmetry
  if (!kramers_coeff)
    scoeff = scoeff->init_kramers_coeff(geom_, overlap_, hcore_, 2*ref_->nclosed() + ref_->nact(), tsymm_, gaunt_, breit_);

  // specify active orbitals and move into the active space
  set<int> active_indices;
  const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
  if (iactive) {
    // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
    for (auto& i : *iactive)
      active_indices.insert(lexical_cast<int>(i->data()) - 1);
    scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_, tsymm_);
  }

  coeff_ = scoeff->block_format();

  mute_stdcout();
  // CASSCF methods should have FCI member. Inserting "ncore" and "norb" keyword for closed and active orbitals.
  if (nact_)
    fci_ = make_shared<ZHarrison>(idata_, geom_, ref_, nclosed_, nact_, coeff_, /*store*/true);
  nstate_ = nact_ ? fci_->nstate() : 1;
  resume_stdcout();

  cout <<  "  === Dirac CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

}


void ZCASSCF::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "      CASSCF calculation     " << endl;
  cout << "  ---------------------------" << endl << endl;
}


void ZCASSCF::print_iteration(const int iter, const vector<double>& energy, const double error, const double time) const {
  if (energy.size() != 1 && iter) cout << endl;
  int i = 0;
  for (auto& e : energy) {
    cout << "     " << setw(5) << iter << setw(4) << i << " " << setw(19) << fixed << setprecision(8) << e << "   "
                    << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2) << time << endl;
    ++i;
  }
}


static streambuf* backup_stream_;
static ofstream* ofs_;

void ZCASSCF::mute_stdcout() const {
  string filename = "casscf.log";
  ofstream* ofs(new ofstream(filename,(backup_stream_ ? ios::app : ios::trunc)));
  ofs_ = ofs;
  backup_stream_ = cout.rdbuf(ofs->rdbuf());
}


void ZCASSCF::resume_stdcout() const {
  cout.rdbuf(backup_stream_);
  delete ofs_;
}


pair<shared_ptr<ZMatrix>, VectorB> ZCASSCF::make_natural_orbitals(shared_ptr<const ZMatrix> rdm1) const {

  // input should be 1rdm in kramers format
  shared_ptr<ZMatrix> tmp;
  if (tsymm_) {
    tmp = make_shared<QuatMatrix>(*rdm1);
#ifndef NDEBUG
    // Failures here can sometimes be fixed by using a tighter convergence threshold in the FCI part
    auto quatrdm = static_pointer_cast<const QuatMatrix>(tmp);
    assert(quatrdm->is_t_symmetric());
#endif
  } else {
    tmp = make_shared<ZMatrix>(*rdm1);
  }

  const bool unitmat = tmp->is_identity(1.0e-14);
  shared_ptr<ZMatrix> out;
  VectorB vec2(tmp->ndim());

  if (!unitmat) {
    VectorB vec(rdm1->ndim());
    tmp->diagonalize(vec);

    if (!tsymm_)
      RelCoeff::rearrange_eig(vec, tmp, false);

    map<int,int> emap;
    out = tmp->clone();

    const bool occ_sort = idata_->get<bool>("occ_sort",false);
    if (tsymm_) {
      if (occ_sort) {
        // sort by natural orbital occupation numbers
        int b2n = out->ndim();
        for (int i = 0; i != out->mdim()/2; ++i) {
          copy_n(tmp->element_ptr(0, out->mdim()/2-1-i), b2n, out->element_ptr(0, i));
          copy_n(tmp->element_ptr(0, out->mdim()-1-i), b2n, out->element_ptr(0, i+b2n/2));
          vec2[b2n/2-i-1] = vec[i] > 0.0 ? vec[i] : 0.0;
          vec2[b2n-1-i] = vec[i] > 0.0 ? vec[i] : 0.0;;
        }
        // fix the phase
        for (int i = 0; i != tmp->ndim(); ++i) {
          if (real(out->element(i,i)) < 0.0)
            blas::scale_n(-1.0, out->element_ptr(0,i), tmp->ndim());
        }
      } else {
        // sort eigenvectors so that buf is close to a unit matrix
        // assumes quaternion symmetry - only the top-left quarter is checked
        // target column
        for (int i = 0; i != tmp->ndim()/2; ++i) {
          // first find the source column
          tuple<int, double> max = make_tuple(-1, 0.0);
          for (int j = 0; j != tmp->ndim()/2; ++j)
            if (std::abs(tmp->element(i,j)) > get<1>(max))
              max = make_tuple(j, std::abs(tmp->element(i,j)));

          // register to emap
          if (emap.find(get<0>(max)) != emap.end()) throw logic_error("In ZCASSCF::make_natural_orbitals(), two columns had max values in the same positions.  This should not happen.");
          assert(get<0>(max) != -1); // can happen if all checked elements are zero, for example
          emap.emplace(get<0>(max), i);

          // copy to the target
          copy_n(tmp->element_ptr(0,get<0>(max)), tmp->ndim(), out->element_ptr(0,i));
          copy_n(tmp->element_ptr(0,get<0>(max)+tmp->ndim()/2), tmp->ndim(), out->element_ptr(0,i+tmp->ndim()/2));
          vec2[i] = vec[get<0>(max)];
          vec2[i+tmp->ndim()/2] = vec[get<0>(max)];
        }

        // fix the phase
        for (int i = 0; i != tmp->ndim(); ++i) {
          if (real(out->element(i,i)) < 0.0)
          blas::scale_n(-1.0, out->element_ptr(0,i), tmp->ndim());
        }
      }
    } else {
      // assumes no particular symmetry - the full matrix is checked
      if (occ_sort) {
        // sort by natural orbital occupation numbers
        int b2n = out->ndim();
        for (int i = 0; i != out->mdim(); ++i) {
          copy_n(tmp->element_ptr(0, out->mdim()-1-i), b2n, out->element_ptr(0, i));
          vec2[b2n-i-1] = vec[i] > 0.0 ? vec[i] : 0.0;
        }
        // fix the phase
        for (int i = 0; i != tmp->ndim(); ++i) {
          const double phase = std::arg(out->element(i,i));
          blas::scale_n(polar(1.0, -phase), out->element_ptr(0,i), tmp->ndim());
        }
      } else {
        // target column
        for (int i = 0; i != tmp->ndim(); ++i) {
          // first find the source column
          tuple<int, double> max = make_tuple(-1, 0.0);
          for (int j = 0; j != tmp->ndim(); ++j)
            if (std::abs(tmp->element(i,j)) > get<1>(max))
              max = make_tuple(j, std::abs(tmp->element(i,j)));

          // register to emap
          if (emap.find(get<0>(max)) != emap.end()) throw logic_error("In ZCASSCF::make_natural_orbitals(), two columns had max values in the same positions.  This should not happen.");
          assert(get<0>(max) != -1); // can happen if all checked elements are zero, for example
          emap.emplace(get<0>(max), i);

          // copy to the target
          copy_n(tmp->element_ptr(0,get<0>(max)), tmp->ndim(), out->element_ptr(0,i));
          vec2[i] = vec[get<0>(max)];
        }

        // fix the phase
        for (int i = 0; i != tmp->ndim(); ++i) {
          const double phase = std::arg(out->element(i,i));
          blas::scale_n(polar(1.0, -phase), out->element_ptr(0,i), tmp->ndim());
        }
      }
    }

  } else { // set occupation numbers, but coefficients don't need to be updated
    for (int i = 0; i != tmp->ndim(); ++i)
      vec2[i] = tmp->get_real_part()->element(i,i);
    out = tmp;
  }

  return make_pair(out, vec2);
}


shared_ptr<const ZMatrix> ZCASSCF::natorb_rdm2_transform(const shared_ptr<ZMatrix> coeff, shared_ptr<const ZMatrix> rdm2) const {
  shared_ptr<ZMatrix> out = rdm2->clone();
  shared_ptr<ZMatrix> buf = rdm2->clone();
  shared_ptr<const ZMatrix> coeff_conjg = coeff->get_conjg();
  const int ndim  = nact_*2;
  const int ndim2 = ndim*ndim;;
  auto half_trans = [&](shared_ptr<const ZMatrix> a, shared_ptr<ZMatrix> b, shared_ptr<ZMatrix> c) {
    zgemm3m_("N", "N", ndim2*ndim, ndim, ndim, 1.0, a->data(), ndim2*ndim, coeff->data(), ndim, 0.0, b->data(), ndim2*ndim);
    for (int i = 0; i != ndim; ++i)
      zgemm3m_("N", "N", ndim2, ndim, ndim, 1.0, b->data()+i*ndim2*ndim, ndim2, coeff_conjg->data(), ndim, 0.0, c->data()+i*ndim2*ndim, ndim2);
  };
  half_trans(rdm2, buf, out);
  blas::transpose(out->data(), ndim2, ndim2, buf->data());
  half_trans(buf, out, buf);
  blas::transpose(buf->data(), ndim2, ndim2, out->data());
  return out;
}


shared_ptr<const RelCoeff_Block> ZCASSCF::update_coeff(shared_ptr<const RelCoeff_Block> cold, shared_ptr<const ZMatrix> natorb) const {
  // D_rs = C*_ri D_ij (C*_rj)^+. Dij = U_ik L_k (U_jk)^+. So, C'_ri = C_ri * U*_ik ; hence conjugation needed
  auto cnew = make_shared<RelCoeff_Block>(*cold, cold->nclosed(), cold->nact(), cold->nvirt_nr(), cold->nneg());
  cnew->copy_block(0, nclosed_*2, cnew->ndim(), nact_*2, cold->slice(nclosed_*2, nocc_*2) * *natorb->get_conjg());
  return cnew;
}


shared_ptr<const ZMatrix> ZCASSCF::update_qvec(shared_ptr<const ZMatrix> qold, shared_ptr<const ZMatrix> natorb) const {
  auto qnew = make_shared<ZMatrix>(*qold * *natorb);
  qnew->copy_block(nclosed_*2, 0, nact_*2, nact_*2, *natorb % *qnew->get_submatrix(nclosed_*2, 0, nact_*2, nact_*2));
  return qnew;
}


shared_ptr<const Reference> ZCASSCF::conv_to_ref() const {
  auto out = make_shared<RelReference>(geom_, coeff_->striped_format(), energy_, nneg_, nclosed_, nact_, nvirt_-nneg_/2, gaunt_, breit_, /*kramers*/true,
                                       fci_->rdm1_av(), fci_->rdm2_av(), fci_->conv_to_ciwfn());
  return out;
}


void ZCASSCF::print_natocc() const {
  assert(occup_.size() > 0);
  cout << "  ========       state-averaged       ======== " << endl;
  cout << "  ======== natural occupation numbers ======== " << endl;
  const int num = tsymm_ ? occup_.size() / 2 : occup_.size();
  for (int i = 0; i != num; ++i)
    cout << setprecision(4) << "   Orbital " << i << " : " << occup_[i] << endl;
  cout << "  ============================================ " << endl;
}


shared_ptr<ZRotFile> ZCASSCF::copy_electronic_rotations(shared_ptr<const ZRotFile> rot) const {
  const int nr_nvirt = nvirt_ - nneg_/2;
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nr_nvirt*2);
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nr_nvirt;   ++j) {
      out->ele_vc(j, i) = rot->ele_vc(j, i);
      out->ele_vc(j + nr_nvirt, i) = rot->ele_vc(j + nvirt_, i);
      out->ele_vc(j, i + nclosed_) = rot->ele_vc(j, i + nclosed_);
      out->ele_vc(j + nr_nvirt, i + nclosed_) = rot->ele_vc(j + nvirt_, i + nclosed_);
    }
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nr_nvirt;   ++j) {
      out->ele_va(j, i) = rot->ele_va(j, i);
      out->ele_va(j + nr_nvirt, i) = rot->ele_va(j + nvirt_, i);
      out->ele_va(j, i + nact_) = rot->ele_va(j, i + nact_);
      out->ele_va(j + nr_nvirt, i + nact_) = rot->ele_va(j + nvirt_, i + nact_);
    }
  for (int i = 0; i != nact_;   ++i)
    for (int j = 0; j != nclosed_; ++j) {
      out->ele_ca(j, i) = rot->ele_ca(j, i);
      out->ele_ca(j + nclosed_, i) = rot->ele_ca(j + nclosed_, i);
      out->ele_ca(j, i + nact_) = rot->ele_ca(j, i + nact_);
      out->ele_ca(j + nclosed_, i + nact_) = rot->ele_ca(j + nclosed_, i + nact_);
    }
  return out;
}


shared_ptr<ZRotFile> ZCASSCF::copy_positronic_rotations(shared_ptr<const ZRotFile> rot) const {
  const int nvirtnr = nvirt_ - nneg_/2;
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nneg_);
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nneg_/2;   ++j) {
      out->ele_vc(j, i) = rot->ele_vc(j + nvirtnr, i);
      out->ele_vc(j, i + nclosed_) = rot->ele_vc(j + nvirtnr, i + nclosed_);
      out->ele_vc(j + nneg_/2, i)  = rot->ele_vc(j + nvirt_ + nvirtnr, i);
      out->ele_vc(j + nneg_/2, i + nclosed_) = rot->ele_vc(j + nvirt_ + nvirtnr, i + nclosed_);
    }
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nneg_/2;   ++j) {
      out->ele_va(j, i) = rot->ele_va(j + nvirtnr, i);
      out->ele_va(j, i + nact_) = rot->ele_va(j + nvirtnr, i + nact_);
      out->ele_va(j + nneg_/2, i) = rot->ele_va(j + nvirt_ + nvirtnr, i);
      out->ele_va(j + nneg_/2, i + nact_) = rot->ele_va(j + nvirt_ + nvirtnr, i + nact_);
    }
  return out;
}


shared_ptr<ZMatrix> ZCASSCF::compute_active_fock(const ZMatView coeff, shared_ptr<const ZMatrix> rdm1) const {
  // calculate S^1/2 of rdm1
  ZMatrix s(*rdm1);
  s.sqrt();
  return make_shared<DFock>(geom_, hcore_->clone(), coeff * *s.get_conjg(), gaunt_, breit_, /*store half*/false, /*robust*/breit_);
}
