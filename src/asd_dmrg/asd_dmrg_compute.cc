//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg_compute.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <memory>

#include <src/util/timer.h>
#include <src/asd_dmrg/asd_dmrg.h>

using namespace std;
using namespace bagel;

void ASD_DMRG::compute() {
  Timer dmrg_timer;

  shared_ptr<DMRG_Block1> left_block, right_block;

  // Seed lattice
  cout << " ===== Start growing DMRG chain =====" << endl;
  {
    shared_ptr<const Reference> ref = multisite_->build_reference(0, vector<bool>(nsites_, true));
    // CI calculation on site 1 with all other sites at meanfield
    left_block = compute_first_block(prepare_growing_input(0), ref);
    left_blocks_.push_back(left_block);
    cout << "  " << print_progress(0, ">>", "..") << setw(16) << dmrg_timer.tick() << endl;
  }

  // Grow lattice
  for (int site = 1; site < nsites_-1; ++site) {
    vector<bool> meanfield(nsites_, true);
    fill_n(meanfield.begin(), site, false);
    shared_ptr<const Reference> ref = multisite_->build_reference(site, meanfield);
    left_block = grow_block(prepare_growing_input(site), ref, left_block, site);
    left_blocks_.push_back(left_block);
    cout << "  " << print_progress(site, ">>", "..") << setw(16) << dmrg_timer.tick() << endl;
  }
  assert(left_blocks_.size() == nsites_-1);

  right_blocks_.resize(nsites_-1);

  cout << endl << " ===== Starting sweeps =====" << endl;

  for (int iter = 0; iter < maxiter_; ++iter) {
  // Start sweeping backwards
    for (int site = nsites_-1; site > 0; --site) {
      left_block = left_blocks_[site-1];
      right_block = (site == nsites_-1) ? nullptr : right_blocks_[nsites_ - site - 2];
      shared_ptr<const Reference> ref = multisite_->build_reference(site, vector<bool>(nsites_, false));

      right_block = decimate_block(prepare_sweeping_input(site), ref, right_block, left_block, site);
      right_blocks_[nsites_ - site - 1] = right_block;
      cout << "  " << print_progress(site, "<<", "<<") << setw(16) << dmrg_timer.tick() << endl;
    }

  // Sweep forwards
    for (int site = 0; site < nsites_-1; ++site) {
      left_block = (site == 0) ? nullptr : left_blocks_[site-1];
      right_block = right_blocks_[nsites_ - site - 2];
      shared_ptr<const Reference> ref = multisite_->build_reference(site, vector<bool>(nsites_, false));

      left_block = decimate_block(prepare_sweeping_input(site), ref, left_block, right_block, site);
      left_blocks_[site] = left_block;
      cout << "  " << print_progress(site, ">>", ">>") << setw(16) << dmrg_timer.tick() << endl;
    }

    bool conv = true;
    bool drop_perturb = true;

    cout << endl;
    for (int i = 0; i < nstate_; ++i) {
      auto mnmx = minmax_element(sweep_energies_[i].begin(), sweep_energies_[i].end());
      const double sweep_average = accumulate(sweep_energies_[i].begin(), sweep_energies_[i].end(), 0.0)/static_cast<double>(sweep_energies_[i].size());
      const double sweep_range = *mnmx.second - *mnmx.first;

      cout << setw(6) << iter << setw(6) << i << setw(22) << setprecision(12) << sweep_average << setw(16) << setprecision(12) << sweep_range
                                                                            << setw(16) << setprecision(12) << energies_[i] - sweep_average << endl;

      if (abs(sweep_range)>thresh_ || energies_[i]-sweep_average>thresh_)
        conv = false;
      if (energies_[i]-sweep_average>perturb_thresh_)
        drop_perturb = false;
      energies_[i] = sweep_average;
      sweep_energies_[i].clear();
    }
    cout << endl;

    if (perturb_!=0.0) {
      conv = false;
      if (drop_perturb) {
        perturb_ *= 0.1;
        if (perturb_<perturb_min_) {
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
}
