//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <cstdio>
#include <src/smith/relcaspt2/RelCASPT2.h>
#include <src/util/math/linearRM.h>
#include <src/prop/pseudospin/pseudospin.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

RelCASPT2::RelCASPT2::RelCASPT2(shared_ptr<const SMITH_Info<std::complex<double>>> ref) : SpinFreeMethod(ref) {
  nstates_ = info_->nact() ? ref->ciwfn()->nstates() : 1;
  auto eig = f1_->diag();
  eig_.resize(eig.size());
  for (int i = 0; i != eig.size(); ++i)
    eig_[i] = real(eig[i]);

  // MS-CASPT2: t2 and s as MultiTensor (t2all, sall)
  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    auto tmp2 = make_shared<MultiTensor>(nstates_);

    for (int j = 0; j != nstates_; ++j)
      if (!info_->sssr() || i == j) {
        (*tmp)[j] = init_amplitude();
        (*tmp2)[j] = init_residual();
      }

    t2all_.push_back(tmp);
    sall_.push_back(tmp2);
    rall_.push_back(tmp2->clone());
  }

  energy_.resize(nstates_);
  pt2energy_.resize(nstates_);
}


void RelCASPT2::RelCASPT2::solve() {
  Timer timer;
  print_iteration();

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  for (int istate = 0; istate != nstates_; ++istate) { //K states
    if (istate > info_->state_begin() || (istate == info_->state_begin() && info_->restart_iter() == 0))
      t2all_[istate]->fac(istate) = 0.0;
    sall_[istate]->fac(istate)  = 0.0;

    for (int jst=0; jst != nstates_; ++jst) { // <jst|
      if (info_->sssr() && jst != istate)
        continue;
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
      while(!sourceq->done())
        sourceq->next_compute();
    }
  }

  // solve linear equation for t amplitudes
  t2all_ = solve_linear(sall_, t2all_);

  // print out energies
  vector<complex<double>> shift_correction(nstates_*nstates_);
  cout << endl;
  for (int istate = 0; istate != nstates_; ++istate) {
    if (info_->shift() == 0.0) {
      pt2energy_[istate] = energy_[istate] + std::real((*eref_)(istate,istate));
      assert(std::abs(std::imag((*eref_)(istate,istate))) < 1.0e-8);
      cout << "    * RelCASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] <<endl;
    } else {
      // will be used in normq
      n = init_residual();
      complex<double> norm = 0.0;
      for (int jstate = 0; jstate != nstates_; ++jstate) {
        if (istate != jstate && info_->shift_diag())
          continue;
        complex<double> nn = 0.0;
        for (int jst = 0; jst != nstates_; ++jst) { // bra
          for (int ist = 0; ist != nstates_; ++ist) { // ket
            if (info_->sssr() && (jst != istate || ist != jstate))
              continue;
            set_rdm(jst, ist);
            t2 = t2all_[jstate]->at(ist);
            shared_ptr<Queue> normq = make_normq(true, jst == ist);
            while (!normq->done())
              normq->next_compute();
            nn += conj(dot_product_transpose(n, t2all_[istate]->at(jst)));
          }
        }
        shift_correction[istate + nstates_*jstate] = nn;
        if (jstate == istate) norm = nn;
      }

      pt2energy_[istate] = energy_[istate] + std::real((*eref_)(istate,istate)) - info_->shift()*std::real(norm);
      assert(std::abs(std::imag((*eref_)(istate,istate))) < 1.0e-8);
      assert(std::abs(std::imag(norm)) < 1.0e-8);
      cout << "    * RelCASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
      cout << "        w/o shift correction  " << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) <<endl;
      cout <<endl;
    }
  }

  // MS-CASPT2
  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<ZMatrix>(nstates_, nstates_);

    if (info_->shift())
      cout << "    MS-RelCASPT2:  Applying levelshift correction to " << (info_->shift_diag() ? "diagonal" : "all" ) <<  " elements of the effective Hamiltonian."  << endl << endl;

    for (int ist = 0; ist != nstates_; ++ist) {
      auto sist = make_shared<MultiTensor>(nstates_);
      for (int jst=0; jst != nstates_; ++jst) {
        if (sall_[ist]->at(jst)) {
          sist->at(jst) = sall_[ist]->at(jst);
        } else {
          set_rdm(jst, ist);
          s = init_residual();
          shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
          while(!sourceq->done())
            sourceq->next_compute();
          sist->at(jst) = s;
        }
      }

      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          (*heff_)(ist, ist) = pt2energy_[ist];
        } else {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > ]
          (*heff_)(jst, ist) = std::conj(dot_product_transpose(sist, t2all_[jst])) + (*eref_)(jst, ist);
          if (!info_->shift_diag())
            (*heff_)(jst, ist) -= shift_correction[jst+nstates_*ist]*info_->shift();
        }
      }
    }
    heff_->hermite();

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-RelCASPT2 Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    heff_->diagonalize(eig);
    copy_n(eig.data(), nstates_, pt2energy_.data());

    // print out the eigen vector
    cout << endl;
    cout << "    * MS-RelCASPT2 rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate)
      cout << "    * MS-RelCASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    cout << endl << endl;
  } else {
    heff_ = make_shared<ZMatrix>(1,1);
    heff_->element(0,0) = 1.0;
  }
  energy_ = pt2energy_;

  // TODO When the Property class is implemented, this should be one
  if (info_->aniso_data()) {
    if (info_->geom()->magnetism()) {
      cout << "  ** Magnetic anisotropy analysis is currently only available for zero-field calculations; sorry." << endl;
    } else {
      shared_ptr<const RelCIWfn> new_ciwfn = info_->do_ms() ? rotate_ciwfn(info_->ciwfn(), *heff_) : info_->ciwfn();
      const int nspin = info_->aniso_data()->get<int>("nspin", nstates_-1);
      Pseudospin ps(nspin, info_->geom(), new_ciwfn, info_->aniso_data());
      ps.compute(energy_, info_->relref()->relcoeff()->block_format()->active_part());
    }
  }
}


// function to solve linear equation
vector<shared_ptr<MultiTensor_<complex<double>>>> RelCASPT2::RelCASPT2::solve_linear(vector<shared_ptr<MultiTensor_<complex<double>>>> s,
                                                                                     vector<shared_ptr<MultiTensor_<complex<double>>>> t) {
  Timer mtimer;
#ifndef DISABLE_SERIALIZATION
  if (info_->restart()) {
    OArchive archive("RelSMITH_info");
    archive << info_;
    mtimer.tick_print("Save RelSMITH info Archive");
  }
#endif

  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;
  for (int i = 0; i != nstates_; ++i) {  // K states
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i] - info_->shift();
    energy_[i] = 0.0;
    // set guess vector
    if (i > info_->state_begin() || (i == info_->state_begin() && info_->restart_iter() == 0))
      t[i]->zero();
    if (s[i]->rms() < 1.0e-15) {
      print_iteration(0, 0.0, 0.0, mtimer.tick());
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      if (i > info_->state_begin() || (i == info_->state_begin() && info_->restart_iter() == 0))
        update_amplitude(t[i], s[i]);
    }

    auto solver = make_shared<LinearRM<MultiTensor, ZMatrix>>(info_->davidson_subspace(), s[i]);
    for (int iter = 0; iter != info_->maxiter(); ++iter) {
      rall_[i]->zero();
      const double norm = t[i]->norm();
      t[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
        for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
          if (info_->sssr() && (jst != i || ist != i))
            continue;
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t[i]->at(ist);
          r = rall_[i]->at(jst);

          // compute residuals named r for each K
          e0_ = e0all_[i] - info_->shift();
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
        }
      }
      // solve using subspace updates
      rall_[i] = solver->compute_residual(t[i], rall_[i]);
      t[i] = solver->civec();
      t2all_[i] = t[i];

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error = rall_[i]->norm() / pow(rall_[i]->size(), 0.25);
      print_iteration(iter, energy_[i], error, mtimer.tick());
      conv = error < info_->thresh();

#ifndef DISABLE_SERIALIZATION
      if (info_->restart() && (conv || (info_->restart_each_iter() && iter > 0))) {
        string arch = "RelCASPT2_";
        if (mpi__->rank() == 0)
          arch += "t2_" + to_string(i) + (conv ? "_converged" : "_iter_" + to_string(iter));
        else
          arch += "temp_trash";
        {
          OArchive archive(arch);
          archive << t2all_[i];
        }
        mpi__->barrier();
        if (conv)
          mtimer.tick_print("Save T-amplitude Archive (RelSMITH)");
        remove("RelCASPT2_temp_trash.archive");
      }
#endif

      // compute delta t2 and update amplitude
      if (!conv) {
        t[i]->zero();
        update_amplitude(t[i], rall_[i]);
      }
      if (conv) break;
    }
    if (i+1 != nstates_) cout << endl;
    converged &= conv;
  }
  print_iteration(!converged);
  return t;
}


void RelCASPT2::RelCASPT2::load_t2all(shared_ptr<MultiTensor> t2in, const int ist) {
  assert(ist <= info_->state_begin());
  t2all_[ist] = t2in;
}


void RelCASPT2::RelCASPT2::solve_gradient(const int targetJ, const int targetI, shared_ptr<const NacmType> nacmtype, const bool nocider) {
  throw std::logic_error("Nuclear gradients not implemented for RelCASPT2");
}


#endif
