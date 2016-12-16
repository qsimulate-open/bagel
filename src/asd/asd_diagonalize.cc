//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_diagonalize.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/asd/asd_base.h>

using namespace std;
using namespace bagel;


void ASD_base::generate_initial_guess(shared_ptr<Matrix> cc, const vector<DimerSubspace_base>& subspaces, int nstates) {
  int trialsize = 0;
  int nguess = nguess_;

  const int subspace_states = accumulate(subspaces.begin(), subspaces.end(), 0,
                                        [] (int x, const DimerSubspace_base& s) { return s.dimerstates()+x; });

  map<double, vector<int>> seeds;
  for(auto& ispace : subspaces) {
    for(int state = ispace.offset(); state < ispace.offset()+ispace.dimerstates(); ++state)
      seeds[denom_[state]].push_back(state);
  }

  while (trialsize < nstates) {
    vector<int> b;
    b.reserve(nguess);

    for (auto& i : seeds) {
      b.insert(b.end(), i.second.begin(), i.second.end());
      if (b.size() >= nguess) break;
    }

    // build matrix
    auto basis = make_shared<Matrix>(dimerstates_, b.size());
    for (int i = 0; i < b.size(); ++i)
      basis->element(b[i], i) = 1.0;

    // build spin operator
    shared_ptr<Matrix> spn = spin_->apply(*basis);
    spn = make_shared<Matrix>( *spn % *basis );
    VectorB spin_values(b.size());
    spn->diagonalize(spin_values);
    const double expected_spin = 0.25 * static_cast<double>(nspin_ * (nspin_ + 2));
    int start, end;
    for (start = 0; start < spin_values.size(); ++start)
      if (fabs(spin_values(start) - expected_spin) < 1.0e-4) break;
    for (end = start; end < spin_values.size(); ++end)
      if (fabs(spin_values(end) - expected_spin) > 1.0e-4) break;

    trialsize = end - start;

    if (trialsize >= nstates) {
      basis = (*basis * *spn).slice_copy(start, end);

      shared_ptr<const Matrix> sigma = apply_hamiltonian(*basis, subspaces_base());
      auto H = make_shared<Matrix>(*sigma % *basis);
      VectorB energies(trialsize);
      H->diagonalize(energies);

      basis = make_shared<Matrix>(*basis * *H);
      for (int i = 0; i < nstates; ++i)
        copy_n(basis->element_ptr(0, i), basis->ndim(), cc->element_ptr(0, i));
    }
    else if (nguess >= subspace_states) {
      throw runtime_error("Requesting more spin allowed states than exist in ASD space");
    }
    else {
      nguess *= 2;
    }
  }
}


shared_ptr<Matrix> ASD_base::apply_hamiltonian(const Matrix& o, const vector<DimerSubspace_base>& subspaces) {
  const int nstates = o.mdim();

  shared_ptr<Matrix> out = o.clone();
  for (auto iAB = subspaces.begin(); iAB != subspaces.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      if (store_matrix_) {
        dgemm_("N", "N", iAB->dimerstates(), nstates, jAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, joff), hamiltonian_->ndim(),
                                                                               o.element_ptr(joff, 0), o.ndim(),
                                                                          1.0, out->element_ptr(ioff, 0), out->ndim());
        dgemm_("T", "N", jAB->dimerstates(), nstates, iAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, joff), hamiltonian_->ndim(),
                                                                               o.element_ptr(ioff, 0), o.ndim(),
                                                                          1.0, out->element_ptr(joff, 0), out->ndim());
      } else {
        shared_ptr<const Matrix> block = couple_blocks(*jAB, *iAB);

        if (block) {
          dgemm_("N", "N", block->ndim(), nstates, block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(joff, 0), o.ndim());
          dgemm_("T", "N", block->mdim(), nstates, block->ndim(), 1.0, block->data(), block->ndim(), o.element_ptr(joff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), o.ndim());
        }
      }
    }

    if (store_matrix_) {
      dgemm_("N", "N", iAB->dimerstates(), nstates, iAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, ioff), hamiltonian_->ndim(),
                                                                             o.element_ptr(ioff, 0), o.ndim(),
                                                                        1.0, out->element_ptr(ioff, 0), out->ndim());
    } else {
      shared_ptr<const Matrix> block = compute_diagonal_block(*iAB);
      dgemm_("N", "N", block->ndim(), nstates, block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), out->ndim());
    }
  }

  return out;
}

// Davidson diagonalize
//  - cc is initial guess on input, eigenvectors on exit
//  - subspaces is the subspaces over which to apply the Hamiltonian
//  - mute is whether to print convergence info
vector<double> ASD_base::diagonalize(shared_ptr<Matrix>& cc, const vector<DimerSubspace_base>& subspaces, const bool mute) {
  Timer asdtime;
  const int nstates = cc->mdim();

  DavidsonDiag<Matrix> davidson(nstates, davidson_subspace_);

  vector<bool> conv(nstates, false);
  vector<double> out(nstates, 0.0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const Matrix> sigma = apply_hamiltonian(*cc, subspaces);

    vector<shared_ptr<const Matrix>> sigman;
    vector<shared_ptr<const Matrix>> ccn;
    for (int i = 0; i != nstates; ++i) {
      if (!conv[i]) {
        sigman.push_back(sigma->slice_copy(i,i+1));
        ccn.push_back(cc->slice_copy(i,i+1));
      }
      else {
        sigman.push_back(shared_ptr<const Matrix>());
        ccn.push_back(shared_ptr<const Matrix>());
      }
    }
    const vector<double> energies = davidson.compute(ccn, sigman);

    // residual
    vector<shared_ptr<Matrix>> errvec = davidson.residual();
    vector<double> errors;
    for (int i = 0; i != nstates; ++i) {
      errors.push_back(errvec.at(i)->rms());
      conv.at(i) = errors.at(i) < thresh_;
    }

    if (any_of(conv.begin(), conv.end(), [] (const bool t) { return (!t); })) {
      for (int ist = 0; ist != nstates; ++ist) {
        if (conv.at(ist)) continue;
        auto tmp_cc = make_shared<Matrix>(dimerstates_, 1);
        double* target_array = tmp_cc->data();
        double* source_array = errvec.at(ist)->data();
        const double en = energies.at(ist);
        for (auto& space : subspaces) {
          for (int i = space.offset(); i != space.offset()+space.dimerstates(); ++i) {
            target_array[i] = source_array[i] / min(en - denom_[i], -0.1);
          }
        }
        spin_->filter(*tmp_cc, nspin_);
        double nrm = tmp_cc->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        tmp_cc->scale(scal);
        copy_n(tmp_cc->data(), dimerstates_, cc->element_ptr(0, ist));
      }
      cc->broadcast();
    }

    if (!mute) {
      if (nstates != 1 && iter) cout << endl;
      for (int i = 0; i != nstates; ++i) {
        cout << setw(7) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                                << setw(17) << fixed << setprecision(8) << energies[i] << "   "
                                << setw(10) << scientific << setprecision(2) << errors[i] << "   "
                                << fixed << setw(10) << setprecision(2) << asdtime.tick() << endl;
      }
    }
    copy(energies.begin(), energies.end(), out.begin());
    if (all_of(conv.begin(), conv.end(), [] (const bool t) { return (t); })) break;
  }

  vector<shared_ptr<Matrix>> eigenstates = davidson.civec();
  for (int i = 0; i != nstates; ++i)
    copy_n(eigenstates.at(i)->data(), dimerstates_, cc->element_ptr(0, i));

  return out;
}
