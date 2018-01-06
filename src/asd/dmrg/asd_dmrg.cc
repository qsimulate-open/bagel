//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_base.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/dmrg/asd_dmrg.h>
#include <src/scf/hf/fock.h>
#include <src/mat1e/overlap.h>

using namespace std;
using namespace bagel;

ASD_DMRG::ASD_DMRG(shared_ptr<const PTree> input, shared_ptr<const Reference> iref) : input_(input) {
  // initialize parameters
  nsites_ = input_->get<int>("nsites", 1);
  if (nsites_ <= 1) throw runtime_error("nsites_ has to be specified and should be no smaller than two");
  nstate_ = input_->get<int>("nstate", 1);
  charge_ = input_->get<int>("charge", 0);
  nspin_ = input_->get<int>("nspin", 0);
  active_sizes_ = input_->get_vector<int>("active_sizes");
  active_electrons_ = input_->get_vector<int>("active_electrons");
  region_sizes_ = input_->get_vector<int>("region_sizes");
  assert(accumulate(region_sizes_.begin(), region_sizes_.end(), 0) == iref->geom()->natom());
  
  ntrunc_ = input_->get<int>("ntrunc", 0);
  if (ntrunc_ == 0) throw runtime_error("ntrunc_ has to be provided");
  thresh_ = input_->get<double>("thresh", 1.0e-6);
  maxiter_ = input_->get<int>("maxiter", 50);
  // perturbation parameters
  perturb_ = input_->get<double>("perturb", 1.0e-3);
  perturb_thresh_ = input_->get<double>("perturb_thresh", 1.0e-4);
  perturb_min_ = input_->get<double>("perturb_min", 1.0e-5);

  down_thresh_ = input_->get<double>("down_thresh", 1.0e-8);
  auto down = input_->get_child_optional("down_sweep_truncs");
  down_sweep_ = static_cast<bool>(down);
  if (down_sweep_)
    down_sweep_truncs_ = input_->get_vector<int>("down_sweep_truncs");

  auto winput = input_->get_child_optional("weights");
  if (winput)
    weights_ = input_->get_vector<double>("weights", nstate_);
  else
    weights_.resize(nstate_, 1.0/static_cast<double>(nstate_));

  energies_.resize(nstate_);
  sweep_energies_.resize(nstate_);

  // reorder coeff to closed-active-virtual
  rearrange_orbitals(iref);
}


void ASD_DMRG::rearrange_orbitals(shared_ptr<const Reference> iref) {
  auto out_coeff = iref->coeff()->clone();

  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  const int nclosed = (iref->geom()->nele() - charge_ - nactele) / 2;
  assert((iref->geom()->nele() - charge_ - nactele) % 2 == 0);
  const int nactive = accumulate(active_sizes_.begin(), active_sizes_.end(), 0);
  const int nvirt = iref->coeff()->mdim() - nclosed - nactive;
  const int nbasis = iref->geom()->nbasis();

  set<int> active_set;
  {
    vector<int> act_vec = input_->get_vector<int>("active_orbitals");
    for_each(act_vec.begin(), act_vec.end(), [](int& x) { --x; });
    active_set.insert(act_vec.begin(), act_vec.end());
  }

  int pos = nclosed;
  for (const int& aorb : active_set)
    copy_n(iref->coeff()->element_ptr(0, aorb), nbasis, out_coeff->element_ptr(0, pos++));

  int closed_position = 0;
  int virt_position = nclosed + nactive;
  for (int iorb = 0; iorb != out_coeff->mdim(); ++iorb)
    if (active_set.count(iorb) == 0)
      copy_n(iref->coeff()->element_ptr(0, iorb), nbasis, out_coeff->element_ptr(0, (closed_position < nclosed ? closed_position++ : virt_position++)));
  assert(closed_position == nclosed);
  assert(virt_position == out_coeff->mdim());
  
  sref_ = make_shared<Reference>(iref->geom(), make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}


void ASD_DMRG::project_active() {

}

/*
  // collect orbital subspaces info
  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  const int nclosed = (iref->geom()->nele() - charge_ - nactele) / 2;
  assert((iref->geom()->nele() - charge_ - nactele) % 2 == 0);
  const int nactive = accumulate(active_sizes_.begin(), active_sizes_.end(), 0);
  const int nvirt = iref->coeff()->mdim() - nclosed - nactive;
  const int nbasis = iref->geom()->nbasis();

  auto active_subspaces = input_->get_child("active_subspaces");
  set<int> active_set;

  auto out_coeff = iref->coeff()->clone();

  if (active_subspaces->size() == nsites_) {
    // active subspaces info is provided, just reorder the orbitals
    
    vector<set<int>> active_orbitals;
    for (auto site : *active_subspaces) {
      vector<int> active_orbs = site->get_vector<int>("");
      for_each(active_orbs.begin(), active_orbs.end(), [](int& x) { --x; });
      active_orbitals.emplace_back(active_orbs.begin(), active_orbs.end());
      active_set.insert(active_orbs.begin(), active_orbs.end());
    }
    assert(active_set.size() == nactive);
    
    auto in_coeff = iref->coeff();
    assert(in_coeff->ndim() == nbasis);
  
    int active_position = nclosed;
    for (auto site : active_orbitals) {
      for (int aorb : site)
        copy_n(in_coeff->element_ptr(0, aorb), nbasis, out_coeff->element_ptr(0, active_position++));
    }
    assert(active_position == nclosed + nactive);
    
  } else if (active_subspaces->size() == 1){
    // if no manually assigned active orbital subspaces, do projection

    for (auto site : *active_subspaces) {
      vector<int> active_vec = site->get_vector<int>("");
      for_each(active_vec.begin(), active_vec.end(), [](int& x) { --x; });
      active_set.insert(active_vec.begin(), active_vec.end());
    }
    assert(active_set.size() == nactive);

    // delocalized active orbital set
    auto act_orbs = make_shared<Matrix>(nbasis, nactive);
    int pos = 0;
    for (auto iorb : active_set)
      copy_n(iref->coeff()->element_ptr(0, iorb), nbasis, act_orbs->element_ptr(0, pos++));
  
    // region bound info
    vector<pair<int, int>> region_bounds;
    int atomstart = 0; int basis_start = 0;
    for (int isize : region_sizes_) {
      int nbasis = 0;
      for (int iatom = atomstart; iatom != atomstart+isize; ++iatom)
        nbasis += iref->geom()->atoms(iatom)->nbasis();
      region_bounds.emplace_back(basis_start, basis_start+nbasis);
      atomstart += isize;
      basis_start += nbasis;
    }
  
    Overlap ovlp(iref->geom());
    Matrix transform(nactive, nactive);
    {
      int ipos = 0;
      for (int isite = 0; isite != nsites_; ++isite) {
        pair<int, int> bounds = region_bounds[isite];
        const int nsub = active_sizes_[isite];
        const int nbasis = bounds.second - bounds.first;
        shared_ptr<const Matrix> sub_ovlp = ovlp.get_submatrix(bounds.first, bounds.first, nbasis, nbasis);
        shared_ptr<const Matrix> sub_actorb = act_orbs->get_submatrix(bounds.first, 0, nbasis, nactive);
        Matrix SVD(*sub_actorb % *sub_ovlp * *sub_actorb);
        VectorB eigs(nactive);
        SVD.diagonalize(eigs);
        const MatView sub_trans = SVD.slice(nactive-nsub, nactive);
        copy_n(sub_trans.data(), sub_trans.size(), transform.element_ptr(0, ipos));
        ipos += nsub;
      }
      assert(ipos == nactive);
      // Lowdin orthonormalization
      Matrix inv_hf(transform % transform);
      inv_hf.inverse_half();
      transform *= inv_hf;
      cout << endl << string(6, '*') << "  If linear dependency is detected, you have to find better active orbitals to do projection  " << string(6, '*') << endl << endl;
    }
    
    // project coeff
    Matrix projected(*act_orbs * transform);
    copy_n(projected.data(), projected.size(), out_coeff->element_ptr(0, nclosed));

  } else {
    throw runtime_error("Must either provide one total active orbital list or assign active orbitals for each fragment");
  }

  int closed_position = 0;
  int virt_position = nclosed + nactive;
  for (int iorb = 0; iorb != out_coeff->mdim(); ++iorb)
    if (active_set.count(iorb) == 0)
      copy_n(iref->coeff()->element_ptr(0, iorb), nbasis, out_coeff->element_ptr(0, (closed_position < nclosed ? closed_position++ : virt_position++)));
  assert(closed_position == nclosed);
  assert(virt_position == out_coeff->mdim());
  
  // canonicalization
  shared_ptr<const Fock<1>> fock;
  {
    shared_ptr<const Matrix> density = iref->coeff()->form_density_rhf(iref->nclosed());
    const MatView ccoeff = iref->coeff()->slice(0, iref->nclosed());
    fock = make_shared<const Fock<1>>(iref->geom(), iref->hcore(), density, ccoeff);
  }
  const MatView active_mos = out_coeff->slice(nclosed, nclosed+nactive);
  Matrix fock_mo(active_mos % *fock * active_mos);
  VectorB eigs(active_mos.mdim());
  shared_ptr<const Matrix> active_transformation = fock_mo.diagonalize_blocks(eigs, active_sizes_);
  const Matrix transformed_mos(active_mos * *active_transformation);
  copy_n(transformed_mos.data(), transformed_mos.size(), out_coeff->element_ptr(0, nclosed));

  sref_ = make_shared<Reference>(iref->geom(), make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}
*/


shared_ptr<Reference> ASD_DMRG::build_reference(const int site, const vector<bool> meanfield) const {
  assert(meanfield.size() == nsites_ && site < nsites_ && site >= 0);

  vector<shared_ptr<const MatView>> closed_orbitals = { make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed())) };
  const int act_start = accumulate(active_sizes_.begin(), active_sizes_.begin()+site, sref_->nclosed());
  const int act_fence = act_start + active_sizes_.at(site);
  const MatView active_orbitals = sref_->coeff()->slice(act_start, act_fence);

  int current = sref_->nclosed();
  for (int i = 0; i != nsites_; ++i) {
    if (meanfield[i] && i != site)
      closed_orbitals.push_back(make_shared<const MatView>(sref_->coeff()->slice(current, current+(active_electrons_.at(i)+1)/2)));
    current += active_sizes_.at(i);
  }

  const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [](int x, shared_ptr<const MatView>m) { return x + m->mdim(); });
  const int nact = active_orbitals.mdim();

  auto out = make_shared<Matrix>(sref_->geom()->nbasis(), nclosed+nact);
  current = 0;
  closed_orbitals.push_back(make_shared<MatView>(active_orbitals));
  for (auto& orbitals : closed_orbitals) {
    copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
    current += orbitals->mdim();
  }

  return make_shared<Reference>(sref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);
}


string ASD_DMRG::print_progress(const int position, const string left_symbol, const string right_symbol) const {
  stringstream out;
  for (int i = 0; i < position; ++i) out << left_symbol << " ";
  out << "** ";
  for (int i = position+1; i < nsites_; ++i) out << right_symbol << " ";

  return out.str();
}


vector<shared_ptr<PTree>> ASD_DMRG::prepare_growing_input(const int site) const {
  vector<shared_ptr<PTree>> out;

  shared_ptr<PTree> base_input = input_->get_child_optional("ras");
  if (!base_input) base_input = make_shared<PTree>();

  base_input->erase("charge");
  base_input->erase("nspin");
  base_input->erase("nstate");

  auto space = input_->get_child_optional("spaces");
  if (!space) throw runtime_error("spaces must be specified");
  else if ( !(space->size()==1 || space->size()==nsites_) )
    throw runtime_error("Must specify either one \"space\" object or one per site");

  auto space_iter = space->begin();
  if (space->size() == nsites_)
    advance(space_iter, site);

  // Quirk of PTree class requires this intermediate step
  auto space_input = make_shared<PTree>(**space_iter);

  for (auto& siter : *space_input) {
    auto tmp = make_shared<PTree>(*base_input);
    tmp->put("charge", siter->get<string>("charge"));
    tmp->put("nspin", siter->get<string>("nspin"));
    tmp->put("nstate", siter->get<string>("nstate"));
    out.push_back(tmp);
  }

  return out;
}


shared_ptr<PTree> ASD_DMRG::prepare_sweeping_input(const int site) const {
  shared_ptr<PTree> out = input_->get_child_optional("ras");
  if (!out) out = make_shared<PTree>();

  out->erase("charge"); out->put("charge", charge_);
  out->erase("nspin");  out->put("nspin", nspin_);
  out->erase("nstate"); out->put("nstate", input_->get<string>("nstate", "1"));

  return out;
}


