//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_compute.cc
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

#include <memory>

#include <src/util/timer.h>
#include <src/asd/dmrg/asd_dmrg.h>

using namespace std;
using namespace bagel;

void ASD_DMRG::sweep() {
  Timer dmrg_timer;

  shared_ptr<DMRG_Block1> left_block, right_block;

  // Seed lattice
  cout << " ===== Start growing DMRG chain =====" << endl;
  {
    shared_ptr<const Reference> ref = build_reference(0, vector<bool>(nsites_, true));
    left_block = compute_first_block(prepare_growing_input(0), ref);
    left_blocks_.push_back(left_block);
    cout << "  " << print_progress(0, ">>", "..") << setw(16) << dmrg_timer.tick() << endl;
  }

  // Grow lattice
  for (int site = 1; site < nsites_-1; ++site) {
    vector<bool> meanfield(nsites_, true);
    fill_n(meanfield.begin(), site, false);
    shared_ptr<const Reference> ref = build_reference(site, meanfield);
    left_block = grow_block(prepare_growing_input(site), ref, left_block, site);
    left_blocks_.push_back(left_block);
    cout << "  " << print_progress(site, ">>", "..") << setw(16) << dmrg_timer.tick() << endl;
  }
  assert(left_blocks_.size() == nsites_-1);

  right_blocks_.resize(nsites_-1);

  cout << endl << " ===== Starting sweeps =====" << endl << endl;

  cout << "  o convergence threshold: " << setw(8) << setprecision(4) << scientific << thresh_ << endl;

  cout << setw(6) << "iter" << setw(6) << "state" << setw(22) << "sweep average" << setw(16) << "sweep range"
                                                                          << setw(16) << "dE average" <<  endl;
  for (int iter = 0; iter < maxiter_; ++iter) {
    // Start sweeping backwards
    for (int site = nsites_-1; site > 0; --site) {
      left_block = left_blocks_[site-1];
      right_block = (site == nsites_-1) ? nullptr : right_blocks_[nsites_ - site - 2];
      shared_ptr<const Reference> ref = build_reference(site, vector<bool>(nsites_, false));

      right_block = decimate_block(prepare_sweeping_input(site), ref, right_block, left_block, site);
      right_blocks_[nsites_ - site - 1] = right_block;
      cout << "  " << print_progress(site, "<<", "<<") << setw(16) << dmrg_timer.tick() << endl;
    }

    // Sweep forwards
    for (int site = 0; site < nsites_-1; ++site) {
      left_block = (site == 0) ? nullptr : left_blocks_[site-1];
      right_block = right_blocks_[nsites_ - site - 2];
      shared_ptr<const Reference> ref = build_reference(site, vector<bool>(nsites_, false));

      left_block = decimate_block(prepare_sweeping_input(site), ref, left_block, right_block, site);
      left_blocks_[site] = left_block;
      cout << "  " << print_progress(site, ">>", ">>") << setw(16) << dmrg_timer.tick() << endl;
    }

    bool conv = (perturb_ < perturb_min_);
    bool drop_perturb = true;

    cout << endl;
    for (int i = 0; i < nstate_; ++i) {
      auto mnmx = minmax_element(sweep_energies_[i].begin(), sweep_energies_[i].end());
      const double sweep_average = accumulate(sweep_energies_[i].begin(), sweep_energies_[i].end(), 0.0)/static_cast<double>(sweep_energies_[i].size());
      const double sweep_range = *mnmx.second - *mnmx.first;

      if (iter != 0)
        cout << setw(6) << iter << setw(6) << i << setw(18) << setprecision(8) << sweep_average << setw(12) << setprecision(8) << sweep_range
                                                                               << setw(12) << setprecision(8) << energies_[i] - sweep_average << endl;
      else
        cout << setw(6) << iter << setw(6) << i << setw(18) << setprecision(8) << sweep_average << setw(12) << setprecision(8) << sweep_range
                                                                               << setw(12) << "---------" << endl;

      conv &= abs(energies_[i] - sweep_average) < thresh_;
      drop_perturb &= abs(energies_[i] - sweep_average) < perturb_thresh_;

      energies_[i] = sweep_average;
      sweep_energies_[i].clear();
    }
    cout << endl;

    if (perturb_ != 0.0) {
      if (drop_perturb) {
        perturb_ *= 0.1;
        if (perturb_ < perturb_min_) {
          cout << "  o perturbation turned off" << endl;
          perturb_ = 0.0;
        }
        else {
          cout << "  o perturbation lowered to " << perturb_ << endl;
        }
      }
    }

    if (conv) {
      cout << "  * Converged!" << endl;
      break;
    }
  }

  // should probably compute some properties before down sweeping

  if (down_sweep_)
    down_sweep();
}

void ASD_DMRG::down_sweep() {
  cout << endl << " ===== Down sweeping =====" << endl;

  remove_if(down_sweep_truncs_.begin(), down_sweep_truncs_.end(), [this] (const int& t) { return t >= ntrunc_; });
  if (!is_sorted(down_sweep_truncs_.rbegin(), down_sweep_truncs_.rend())) {
    cout << "  o Sorting list of truncations into descending order. Was there some reason to have them unordered?" << endl;
    sort(down_sweep_truncs_.rbegin(), down_sweep_truncs_.rend());
  }

  Timer dmrg_timer(0);

  vector<tuple<int, vector<double>>> trunc_convergence;
  trunc_convergence.emplace_back(ntrunc_, energies_);

  perturb_ = 0.0;

  for (int ntrunc : down_sweep_truncs_) {
    ntrunc_ = ntrunc;

    vector<double> energies(nstate_, 0.0);

    shared_ptr<DMRG_Block1> left_block, right_block;

    cout << "  o Starting sweep with M = " << ntrunc_ << endl;
    for (int iter = 0; iter < maxiter_; ++iter) {
    // Start sweeping backwards
      for (int site = nsites_-1; site > 0; --site) {
        left_block = left_blocks_[site-1];
        right_block = (site == nsites_-1) ? nullptr : right_blocks_[nsites_ - site - 2];
        shared_ptr<const Reference> ref = build_reference(site, vector<bool>(nsites_, false));

        right_block = decimate_block(prepare_sweeping_input(site), ref, right_block, left_block, site);
        right_blocks_[nsites_ - site - 1] = right_block;
        cout << "  " << print_progress(site, "<<", "<<") << setw(16) << dmrg_timer.tick() << endl;
      }

    // Sweep forwards
      for (int site = 0; site < nsites_-1; ++site) {
        left_block = (site == 0) ? nullptr : left_blocks_[site-1];
        right_block = right_blocks_[nsites_ - site - 2];
        shared_ptr<const Reference> ref = build_reference(site, vector<bool>(nsites_, false));

        left_block = decimate_block(prepare_sweeping_input(site), ref, left_block, right_block, site);
        left_blocks_[site] = left_block;
        cout << "  " << print_progress(site, ">>", ">>") << setw(16) << dmrg_timer.tick() << endl;
      }

      bool conv = true;
      cout << endl;
      for (int i = 0; i < nstate_; ++i) {
        auto mnmx = minmax_element(sweep_energies_[i].begin(), sweep_energies_[i].end());
        const double sweep_average = accumulate(sweep_energies_[i].begin(), sweep_energies_[i].end(), 0.0)/static_cast<double>(sweep_energies_[i].size());
        const double sweep_range = *mnmx.second - *mnmx.first;

        if (iter != 0)
          cout << setw(6) << iter << setw(6) << i << setw(18) << setprecision(8) << sweep_average << setw(12) << setprecision(8) << sweep_range
                                                                                  << setw(12) << setprecision(8) << energies[i] - sweep_average << endl;
        else
          cout << setw(6) << iter << setw(6) << i << setw(18) << setprecision(8) << sweep_average << setw(12) << setprecision(8) << sweep_range
                                                                                  << setw(12) << "---------" << endl;

        conv &= abs(energies[i]-sweep_average) < down_thresh_;

        energies[i] = sweep_average;
        sweep_energies_[i].clear();
      }
      cout << endl;

      if (conv) {
        cout << "  * Converged!" << endl;
        break;
      }
    }

    trunc_convergence.emplace_back(ntrunc_, energies);
  }

  cout << endl;

  cout << " Convergence with M:" << endl;

  cout << setw(6) << "M";
  for (int ist = 0; ist < nstate_; ++ist)
    cout << setw(22) << string("state ") + to_string(ist);
  cout << endl;

  for (auto& i : trunc_convergence) {
    cout << setw(6) << get<0>(i);
    for (auto& j : get<1>(i))
      cout << setw(22) << setprecision(16) << j;
    cout << endl;
  }
}
