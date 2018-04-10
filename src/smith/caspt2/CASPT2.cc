//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
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


#include <iostream>
#include <iomanip>
#include <src/smith/caspt2/CASPT2.h>
#include <src/util/math/linearRM.h>
#include <src/smith/caspt2/MSCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = info_->nact() ? ref->ciwfn()->nstates() : 1;

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


CASPT2::CASPT2::CASPT2(const CASPT2& cas) : SpinFreeMethod(cas) {
  info_    = cas.info_;
  virt_    = cas.virt_;
  active_  = cas.active_;
  closed_  = cas.closed_;
  rvirt_   = cas.rvirt_;
  ractive_ = cas.ractive_;
  rclosed_ = cas.rclosed_;
  nstates_ = cas.nstates_;
  heff_    = cas.heff_;
  fockact_ = cas.fockact_;
  e0all_   = cas.e0all_;
  xmsmat_  = cas.xmsmat_;
  energy_  = cas.energy_;
  pt2energy_ = cas.pt2energy_;

  // sall is changed in gradient and nacme codes while the others are not
  t2all_ = cas.t2all_;
  t2all_orthogonal_ = cas.t2all_orthogonal_;
  rall_  = cas.rall_;
  for (int i = 0; i != nstates_; ++i) {
    sall_.push_back(cas.sall_[i]->copy());
  }
  h1_ = cas.h1_;
  f1_ = cas.f1_;
  v2_ = cas.v2_;

  rdm0all_ = cas.rdm0all_;
  rdm1all_ = cas.rdm1all_;
  rdm2all_ = cas.rdm2all_;
  rdm3all_ = cas.rdm3all_;
  rdm4all_ = cas.rdm4all_;
}


void CASPT2::CASPT2::do_rdm_deriv(double factor) {
  Timer timer(1);
  tie(den0cirdmt, den1cirdmt, den2cirdmt, den3cirdmt, den4cirdmt) = feed_denci();

  ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), 1);
  const size_t nact  = info_->nact();
  const size_t norb2 = nact * nact;
  const size_t ndet = ci_deriv_->data(0)->size();
  const size_t ijmax = info_->cimaxchunk();
  const size_t ijnum = ndet * norb2 * norb2;
  const size_t npass = ((mpi__->size() > ((ijnum - 1)/ijmax + 1)) && (mpi__->size() != 1) && ndet > 10000) ? mpi__->size() : (ijnum - 1) / ijmax + 1;
  const size_t nsize = (ndet - 1) / npass + 1;

  if (npass > 1)
    cout << "       - CI derivative contraction will be done with " << npass << " passes" << endl;

  // Fock-weighted 2RDM derivative evaluated first (needed for calculating Fock-weighted 3RDM derivative)
  rdm2fderiv_ = SpinFreeMethod<double>::feed_rdm_2fderiv(info_, fockact_, 0);

  if (npass > 1)
    timer.tick_print("Fock-weighted 2RDM derivative");

  // embarrasingly parallel mode. npass > 1 -> distribute among the nodes.
  // otherwise just do using all the nodes.
  const int nproc = npass > 1 ? 1 : mpi__->size();
  const int ncomm = mpi__->size() / nproc;
  const int icomm = mpi__->rank() / nproc;
  mpi__->split(nproc);

  for (int ipass = 0; ipass != npass; ++ipass) {
    if (ipass % ncomm == icomm && ncomm != icomm) {
      const size_t ioffset = ipass * nsize;
      const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;
      tie(rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3fderiv_)
        = SpinFreeMethod<double>::feed_rdm_deriv(info_, fockact_, 0, ioffset, isize, rdm2fderiv_);

      shared_ptr<VectorB> bdata = contract_rdm_deriv(info_->ciwfn(), ioffset, isize, fockact_);
      blas::ax_plus_y_n(factor, bdata->data(), ndet, ci_deriv_->data(0)->data());

      if (npass > 1) {
        stringstream ss; ss << "Multipassing (" << setw(2) << ipass + 1 << " / " << npass << ")";
        timer.tick_print(ss.str());
      }
    }
  }
  mpi__->merge();

  if (npass > 1)
    mpi__->allreduce(ci_deriv_->data(0)->data(), ndet);
}


void CASPT2::CASPT2::solve() {
  if (info_->orthogonal_basis())
    cout << "    * CASPT2 iteration is performed using orthogonal basis" << endl << endl;
  else {
    cout << "    * CASPT2 iteration is performed using redundant basis" << endl << endl;
  }
  Timer timer;
  if (!info_->orthogonal_basis())
    print_iteration();

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  for (int istate = 0; istate != nstates_; ++istate) { //K states
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
  if (info_->orthogonal_basis()) {
    tie(t2all_orthogonal_, t2all_) = solve_linear_orthogonal(sall_, t2all_);
  } else {
    t2all_ = solve_linear(sall_, t2all_);
  }

  timer.tick_print("CASPT2 energy evaluation");
  cout << endl;

  // recalculate the energy without the shift
  for (int istate = 0; istate != nstates_; ++istate) {
    rall_[istate]->zero();
    for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
      for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
        if (info_->sssr() && (jst != istate || ist != istate))
          continue;
        set_rdm(jst, ist);
        t2 = t2all_[istate]->at(ist);
        r = rall_[istate]->at(jst);

        e0_ = e0all_[istate];
        shared_ptr<Queue> queue = make_residualq(true, jst == ist);
        while (!queue->done())
          queue->next_compute();
        diagonal(r, t2, jst == ist);
      }
    }
    pt2energy_[istate] = (*eref_)(istate,istate) + detail::real(dot_product_transpose(sall_[istate], t2all_[istate])) * 2.0;
    pt2energy_[istate] += detail::real(dot_product_transpose(rall_[istate], t2all_[istate]));
    cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    if (info_->shift() != 0.0)
      cout << "        w/o shift correction  " << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) << endl;
    {
      n = init_residual();
      double norm = 0.0;
      for (int jst = 0; jst != nstates_; ++jst) { // bra
        for (int ist = 0; ist != nstates_; ++ist) { // ket
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          shared_ptr<Queue> normq = make_normq(true, jst == ist);
          while (!normq->done())
            normq->next_compute();
          norm += dot_product_transpose(n, t2all_[istate]->at(jst));
        }
      }
      cout << setprecision(10) << "norm = " << norm << endl;
      cout << "        reference weight      " << fixed << setw(20) << setprecision(10) << 1.0/(1.0+norm) << endl;
    }
    cout << endl;
  }

  // TODO Implement off-diagonal shift correction for nonrelativistic energy + nuclear gradients
  if (info_->shift() && info_->do_ms() && !info_->shift_diag())
    cout << "    Applying levelshift correction to diagonal elements of the Hamiltonian only.  (Off-diagonals have not been implemented for non-relativistic CASPT2.)" << endl << endl;

  // MS-CASPT2
  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<Matrix>(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      auto sist = make_shared<MultiTensor>(nstates_);
      for (int jst = 0; jst != nstates_; ++jst) {
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
          (*heff_)(jst, ist) = dot_product_transpose(sist, t2all_[jst]) + (*eref_)(jst, ist);
        }
      }
    }
    heff_->symmetrize();

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASPT2 Heff";
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
    cout << "    * MS-CASPT2 rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    if (xmsmat_) {
      cout << endl;
      cout << "    * XMS-CASPT2 rotation matrix";
      for (int ist = 0; ist != nstates_; ++ist) {
        cout << endl << "      ";
        for (int jst = 0; jst != nstates_; ++jst)
          cout << setw(20) << setprecision(10) << msrot()->element(ist, jst);
      }
      cout << endl << endl;
    }

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate)
      cout << "    * MS-CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    cout << endl << endl;
  } else {
    heff_ = make_shared<Matrix>(1,1);
    heff_->element(0,0) = 1.0;
  }
  energy_ = pt2energy_;
}


// temporary
void CASPT2::CASPT2::manipulate(shared_ptr<MultiTensor_<double>> s) {
  for (int i = 0; i != nstates_; ++i) {
    if (!s->at(i)) continue;
    // a i b j
#if 0
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = s->at(i)->get_size(i0, i1, i2, i3);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i0, i1, i2, i3);
          }
#endif
    // a r b s
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i3 : virt_)
          for (auto& i1 : virt_) {
            const size_t blocksize = s->at(i)->get_size(i0, i1, i2, i3);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i0, i1, i2, i3);
          }
    // a i r j
    for (auto& i3 : active_) 
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = s->at(i)->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i2, i3, i0, i1);
          }
    // a r b i
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i2 : closed_)
          for (auto& i1 : virt_) {
            const size_t blocksize = s->at(i)->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i2, i3, i0, i1);
          }
    // r i s j
    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i2 : closed_)
          for (auto& i0 : closed_) {
            const size_t blocksize = s->at(i)->get_size(i0, i1, i2, i3);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i0, i1, i2, i3);
          }
    // a i r s & a r s i
    for (auto& i3 : active_)
      for (auto& i2 : active_)
        for (auto& i1 : virt_)
          for (auto& i0 : closed_) {
            const size_t blocksize = s->at(i)->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i2, i3, i0, i1);
            s->at(i)->put_block(data, i0, i3, i2, i1);
          }
    // a r s t
    for (auto& i3 : active_)
      for (auto& i2 : active_)
        for (auto& i0 : active_)
          for (auto& i1 : virt_) {
            const size_t blocksize = s->at(i)->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i2, i3, i0, i1);
          }
    // r i s t
    for (auto& i3 : active_)
      for (auto& i1 : active_)
        for (auto& i0 : active_)
          for (auto& i2 : closed_) {
            const size_t blocksize = s->at(i)->get_size(i2, i3, i0, i1);
            unique_ptr<double[]> data(new double[blocksize]);
            for (int j = 0; j != blocksize; ++j) data[j] = 0.0;
            s->at(i)->put_block(data, i2, i3, i0, i1);
          }
  }
}

// function to solve linear equation
vector<shared_ptr<MultiTensor_<double>>> CASPT2::CASPT2::solve_linear(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  Timer mtimer;
  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;
  for (int i = 0; i != nstates_; ++i) {  // K states
    ///vvvvvvv
    manipulate(s[i]);
    ///^^^^^^^
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i] - info_->shift();
    energy_[i] = 0.0;
    // set guess vector
    t[i]->zero();
    if (s[i]->rms() < 1.0e-15) {
      print_iteration(0, 0.0, 0.0, mtimer.tick());
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      update_amplitude(t[i], s[i]);
    }

    auto solver = make_shared<LinearRM<MultiTensor>>(info_->davidson_subspace(), s[i]);
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
      ///vvvvvvv
      manipulate(rall_[i]);
      ///^^^^^^^

      // solve using subspace updates
      rall_[i] = solver->compute_residual(t[i], rall_[i]);
      t[i] = solver->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error = rall_[i]->norm() / pow(rall_[i]->size(), 0.25);
      print_iteration(iter, energy_[i], error, mtimer.tick());
      conv = error < info_->thresh();

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


tuple<vector<shared_ptr<VectorB>>,vector<shared_ptr<MultiTensor_<double>>>>
CASPT2::CASPT2::solve_linear_orthogonal(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  Timer mtimer;
  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;
  vector<shared_ptr<VectorB>> out;
  cout << endl << "      -----------------------------------------------------  CASPT2 iteration  --------------------------------------------------------------" << endl;
  cout << "       #        aibj        arbs        arbi        airj        risj        airs        arst        rist           Etot          error   time" << endl;
  cout << "      ---------------------------------------------------------------------------------------------------------------------------------------" << endl << endl;
  for (int i = 0; i != nstates_; ++i) {  // K states
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i];
    energy_[i] = 0.0;
    ///vvvvvvv
    manipulate(s[i]);
    ///^^^^^^^
    shared_ptr<VectorB> source = transform_to_orthogonal(s[i]);
    // set guess vector
    auto amplitude = make_shared<VectorB>(source->size());
    if (s[i]->rms() < 1.0e-15) {
      print_energy_parts(0, source, source, amplitude, 0.0, mtimer.tick());
      out.push_back(amplitude);
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      update_amplitude_orthogonal(amplitude, source, nstates_, i);
    }

    // TODO should replace (or distribute) the parameters among the nodes
    auto solver = make_shared<LinearRM<VectorB>>(info_->davidson_subspace(), source);
    for (int iter = 0; iter != info_->maxiter(); ++iter) {
      rall_[i]->zero();

      const double norm = amplitude->norm();
      amplitude->scale(1.0/norm);
      t[i] = transform_to_redundant_amplitude(amplitude, nstates_, i);

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
          e0_ = e0all_[i];
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();

          diagonal(r, t2, jst == ist);
        }
      }
      ///vvvvvvv
      manipulate(rall_[i]);
      ///^^^^^^^
      shared_ptr<VectorB> residual = transform_to_orthogonal(rall_[i]);
      if (info_->shift() != 0.0)
        add_shift(residual, amplitude, i);
      residual = solver->compute_residual(amplitude, residual);
      amplitude = solver->civec();

      // compute rms for state i
      error = residual->norm() / pow(residual->size(), 0.25);
      energy_[i] = print_energy_parts(iter, source, residual, amplitude, error, mtimer.tick());
      conv = error < info_->thresh();

      if (!conv) {
        amplitude->zero();
        update_amplitude_orthogonal(amplitude, residual, nstates_, i);
      }
      if (conv) {
        t[i] = transform_to_redundant_amplitude(amplitude, nstates_, i);
        const double ntnt = compute_norm(amplitude, amplitude);
        cout << setprecision(10) << "norm = " << ntnt << endl;
        out.push_back(amplitude);
        break;
      }
    }
    if (i+1 != nstates_) cout << endl;
    converged &= conv;
  }
  cout << endl << "      ---------------------------------------------------------------------------------------------------------------------------------------" << endl << endl;
  return make_tuple(out, t);
}


void CASPT2::CASPT2::solve_dm(const int istate, const int jstate) {
  {
    MSCASPT2::MSCASPT2 ms(*this);
    ms.solve_dm(istate, jstate);
    vden1_ = ms.vden1();
  }
}


void CASPT2::CASPT2::solve_gradient(const int targetJ, const int targetI, shared_ptr<const NacmType> nacmtype, const bool nocider) {
  Timer timer;
  // First solve lambda equation if this is MS-CASPT2
  assert (!((targetJ != targetI) && (nstates_ == 1)));

  if ((info_->do_ms() && nstates_ > 1) || info_->shift() != 0.0) {
    // Lambda equation solver
    for (int i = 0; i != nstates_; ++i)
      lall_.push_back(t2all_[i]->clone());
    if (!info_->orthogonal_basis())
      print_iteration();

    // source stores the result of summation over M'
    if (targetJ == targetI) {
      // Gradient: is special case of targetJ = targetI.
      const int target = targetJ;

      auto source = make_shared<MultiTensor>(nstates_);
      for (auto& i : *source)
        i = init_residual();
      if (info_->shift_imag()) {
        // This should also yield the right results for the real shift.
        // Nevertheless, it requires additional evaluation of residual-like term,
        // and therefore, only applied for imaginary case only
        for (int ist = 0; ist != nstates_; ++ist) { // N states
          sall_[ist]->zero();
          auto sist = make_shared<MultiTensor>(nstates_);
          for (int jst = 0; jst != nstates_; ++jst) { // M states
            set_rdm(jst, ist);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
            while(!sourceq->done())
              sourceq->next_compute();
            sist->at(jst) = s;
          }
          source->ax_plus_y((*heff_)(ist, target), sist);
        }

        for (int istate = 0; istate != nstates_; ++istate) { //L states
          for (int jst = 0; jst != nstates_; ++jst) { // M
            if (!info_->sssr() || istate == jst)
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, target), source->at(jst));
          }
        }

        for (int istate = 0; istate != nstates_; ++istate) {  // L states
          for (int jst = 0; jst != nstates_; ++jst) {
            if (info_->sssr() && istate != jst) continue;
            set_rdm(jst, istate);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
            while(!sourceq->done())
              sourceq->next_compute();
            sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, target) * (*heff_)(istate, target), s);
          }
        }

        for (int i = 0; i != nstates_; ++i) {             // L states
          for (int jst = 0; jst != nstates_; ++jst) {   // M states
            for (int ist = 0; ist != nstates_; ++ist) {     // N states
              if (info_->sssr() && (jst != i || ist != i))
                continue;
              set_rdm(jst, ist);
              t2 = t2all_[i]->at(ist);
              r = init_residual();
              e0_ = e0all_[i];
              shared_ptr<Queue> queue = make_residualq(true, jst == ist);
              while (!queue->done())
                queue->next_compute();
              diagonal(r, t2, ist == jst);

              sall_[i]->at(jst)->ax_plus_y(2.0 * (*heff_)(i, target) * (*heff_)(i, target), r);
            }
          }
        }
      } else {
        for (int ist = 0; ist != nstates_; ++ist) {//N states
          auto sist = make_shared<MultiTensor>(nstates_);
          for (int jst = 0; jst != nstates_; ++jst) {
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
          source->ax_plus_y((*heff_)(ist, target), sist);
        }

        for (int istate = 0; istate != nstates_; ++istate) { //L states
          sall_[istate]->zero();
          for (int jst = 0; jst != nstates_; ++jst)
            if (!info_->sssr() || istate == jst)
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, target), source->at(jst));
          if (info_->shift() != 0.0) {
            // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
            n = init_residual();
            for (int jst = 0; jst != nstates_; ++jst) { // bra
              for (int ist = 0; ist != nstates_; ++ist) { // ket
                if (info_->sssr() && (jst != istate || ist != istate))
                  continue;
                set_rdm(jst, ist);
                t2 = t2all_[istate]->at(ist);
                shared_ptr<Queue> normq = make_normq(true, jst == ist);
                while (!normq->done())
                  normq->next_compute();
                sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * pow((*heff_)(istate, target), 2.0), n);
              }
            }
          }
        }
      }
    } else {
      auto sourceJ = make_shared<MultiTensor>(nstates_);
      auto sourceI = make_shared<MultiTensor>(nstates_);
      for (auto& i : *sourceJ)
        i = init_residual();
      for (auto& i : *sourceI)
        i = init_residual();
      // NACME case
      if (info_->shift_imag()) {
        // This should also yield the right results for the real shift.
        // Nevertheless, it requires additional evaluation of residual-like term,
        // and therefore, only applied for imaginary case only
        for (int ist = 0; ist != nstates_; ++ist) { // L states
          sall_[ist]->zero();
          auto sist = make_shared<MultiTensor>(nstates_);
          for (int jst = 0; jst != nstates_; ++jst) {
            set_rdm(jst, ist);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
            while(!sourceq->done())
              sourceq->next_compute();
            sist->at(jst) = s;
          }
          sourceJ->ax_plus_y((*heff_)(ist, targetI) * 0.5, sist);
          sourceI->ax_plus_y((*heff_)(ist, targetJ) * 0.5, sist);
        }

        for (int istate = 0; istate != nstates_; ++istate) { //L states
          for (int jst = 0; jst != nstates_; ++jst) { // M states
            if (!info_->sssr() || istate == jst) {
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetI), sourceI->at(jst));
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetJ), sourceJ->at(jst));
            }
          }
        }

        for (int istate = 0; istate != nstates_; ++istate) {  // L states
          for (int jst = 0; jst != nstates_; ++jst) {
            if (info_->sssr() && istate != jst) continue;
            set_rdm(jst, istate);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
            while(!sourceq->done())
              sourceq->next_compute();
            sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetI) * (*heff_)(istate, targetJ), s);
          }
        }

        for (int i = 0; i != nstates_; ++i) {             // L states
          for (int jst = 0; jst != nstates_; ++jst) {     // N states
            for (int ist = 0; ist != nstates_; ++ist) {   // M states
              if (info_->sssr() && (jst != i || ist != i))
                continue;
              set_rdm(jst, ist);
              t2 = t2all_[i]->at(ist);
              r = init_residual();
              e0_ = e0all_[i];
              shared_ptr<Queue> queue = make_residualq(true, ist == jst);
              while (!queue->done())
                queue->next_compute();
              diagonal(r, t2, ist == jst);

              sall_[i]->at(jst)->ax_plus_y(2.0 * (*heff_)(i, targetJ) * (*heff_)(i, targetI), r);
            }
          }
        }
      } else {
        for (int ist = 0; ist != nstates_; ++ist) { // L states
          auto sist = make_shared<MultiTensor>(nstates_);
          for (int jst = 0; jst != nstates_; ++jst) {
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
          sourceJ->ax_plus_y((*heff_)(ist, targetI) * 0.5, sist);
          sourceI->ax_plus_y((*heff_)(ist, targetJ) * 0.5, sist);
        }

        for (int istate = 0; istate != nstates_; ++istate) { //K states
          sall_[istate]->zero();
          for (int jst = 0; jst != nstates_; ++jst)
            if (!info_->sssr() || istate == jst) {
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetI), sourceI->at(jst));
              sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetJ), sourceJ->at(jst));
            }
          if (info_->shift() != 0.0) {
            // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
            n = init_residual();
            for (int jst = 0; jst != nstates_; ++jst) { // bra
              for (int ist = 0; ist != nstates_; ++ist) { // ket
                if (info_->sssr() && (jst != istate || ist != istate))
                  continue;
                set_rdm(jst, ist);
                t2 = t2all_[istate]->at(ist);
                shared_ptr<Queue> normq = make_normq(true, jst == ist);
                while (!normq->done())
                  normq->next_compute();
                sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * (*heff_)(istate, targetJ) * (*heff_)(istate, targetI), n);
              }
            }
          }
        }
      }
    }
    // solve linear equation and store lambda in lall
    if (info_->orthogonal_basis()) {
      tie(lall_orthogonal_, lall_) = solve_linear_orthogonal(sall_, lall_);
    } else {
      lall_ = solve_linear(sall_, lall_);
    }

    timer.tick_print("CASPT2 lambda equation");
  }

  if (lall_.empty()) {
    t2 = t2all_[0]->at(0);
    {
      den2 = h1_->clone();
      shared_ptr<Queue> dens2 = make_densityq();
      while (!dens2->done())
        dens2->next_compute();
      den2_ = den2->matrix();
    } {
      den1 = h1_->clone();
      shared_ptr<Queue> dens1 = make_density1q();
      while (!dens1->done())
        dens1->next_compute();
      den1_ = den1->matrix();
    } {
      Den1 = init_residual();
      shared_ptr<Queue> Dens1 = make_density2q();
      while (!Dens1->done())
        Dens1->next_compute();
      Den1_ = Den1;
    }
    timer.tick_print("Correlated density matrix evaluation");

    // then form deci0 - 4
    den0ci = rdm0_->clone();
    den1ci = rdm1_->clone();
    den2ci = rdm2_->clone();
    den3ci = rdm3_->clone();
    den4ci = rdm3_->clone();
    shared_ptr<Queue> dec = make_deciq(/*reset = */true);
    while (!dec->done())
      dec->next_compute();
    timer.tick_print("CI derivative evaluation");

    // when active is divided into the blocks, den4ci is evaluated (activeblock)**2 times
    double den4factor = 1.0 / static_cast<double>(active_.nblock() * active_.nblock());
    den4ci->scale(den4factor);

    // and contract them with rdm derivs
    do_rdm_deriv(1.0);

    timer.tick_print("CI derivative contraction");
    cout << endl;
  } else {
    // in case when CASPT2 is not variational...
    MSCASPT2::MSCASPT2 ms(*this);
    ms.solve_gradient(targetJ, targetI, nocider);
    den1_ = ms.rdm11();
    den2_ = ms.rdm12();
    if (info_->shift_imag())
      den2_tt_ = ms.rdm12_tt();
    Den1_ = ms.rdm21();
    if (!nocider)
      ci_deriv_ = ms.ci_deriv();
    dcheck_ = ms.dcheck();
    if (targetJ != targetI)
      vden1_ = ms.vden1();
    timer.tick();
  }

  correlated_norm_lt_.resize(nstates_);
  correlated_norm_tt_.resize(nstates_);
  if (nstates_ == 1 && info_->shift() == 0.0) {
    n = init_residual();
    shared_ptr<Queue> normq = make_normq();
    while (!normq->done())
      normq->next_compute();
    correlated_norm_lt_[0] = dot_product_transpose(n, t2);
  } else {
    n = init_residual();
    for (int istate = 0; istate != nstates_; ++istate) {
      double tmp = 0.0;
      double tmp2 = 0.0;
      for (int jst = 0; jst != nstates_; ++jst) { // bra
        for (int ist = 0; ist != nstates_; ++ist) { // ket
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          shared_ptr<Queue> normq = make_normq(true, jst == ist);
          while (!normq->done())
            normq->next_compute();
          tmp += dot_product_transpose(n, lall_[istate]->at(jst));
          tmp2 += dot_product_transpose(n, t2all_[istate]->at(jst));
        }
      }
      correlated_norm_lt_[istate] = tmp;
      correlated_norm_tt_[istate] = tmp2;
    }
  }
  timer.tick_print("T1 norm evaluation");

  // some additional terms to be added
  const int ncore = info_->ncore();
  const int nclosed = info_->nclosed()-info_->ncore();
  const int nact = info_->nact();

  // TODO cleanup
#if 1
  if (info_->shift_imag()) {
    const size_t nact = info_->nact();
    const size_t nclosed = info_->nclosed();
    const size_t nvirt = info_->nvirt();
    const size_t nocc = nact + nclosed;
    const size_t ncore = info_->ncore();
    const size_t nclo = nclosed - ncore;
    const size_t size_aibj = nvirt * nvirt * nclo * nclo;
    const size_t size_arbs = denom_->shalf_xx()->ndim()  * nvirt * nvirt;
    const size_t size_arbi = denom_->shalf_x()->ndim()   * nvirt * nclo * nvirt;
    const size_t size_airj = denom_->shalf_h()->ndim()   * nclo * nvirt * nclo;
    const size_t size_risj = denom_->shalf_hh()->ndim()  * nclo * nclo;
    const size_t size_airs = denom_->shalf_xh()->ndim()  * nclo * nvirt;
    const size_t size_arst = denom_->shalf_xxh()->ndim() * nvirt;
    const size_t size_rist = denom_->shalf_xhh()->ndim() * nclo;

    const size_t size_all = size_aibj + size_arbs + size_arbi + size_airj + size_risj + size_airs + size_arst + size_rist;
    const double shift2 = info_->shift() * info_->shift();
    shared_ptr<Matrix> dshift = den2_->clone();
    for (int istate = 0; istate != nstates_; ++istate) { // state of T
      // considering only SS-SR case
      shared_ptr<VectorB> lambda = lall_orthogonal_[istate];
      shared_ptr<VectorB> amplitude = t2all_orthogonal_[istate];
      // a i b j
      {
        size_t ioffset = 0;
        for (int j3 = 0; j3 != nvirt; ++j3)
          for (int j2 = 0; j2 != nclo; ++j2)
            for (int j1 = 0; j1 != nvirt; ++j1)
              for (int j0 = 0; j0 != nclo; ++j0) {
                const int j0i = j0;
                const int j1i = j1 + nocc - ncore;
                const int j2i = j2;
                const int j3i = j3 + nocc - ncore;
                const size_t jall = j0 + nclo * (j1 + nvirt * (j2 + nclo * j3)) + ioffset;
                const size_t jall2 = j0 + nclo * (j3 + nvirt * (j2 + nclo * j1)) + ioffset;
                const double lcovar = ((*lambda)[jall] * 2.0 - (*lambda)[jall2]);
                const double denom = - eig_[j0+ncore] - eig_[j2+ncore] + eig_[j3+nocc] + eig_[j1+nocc];
                // almost all contributions are from closed
                dshift->element(j0i, j0i) += (shift2) * lcovar * (*amplitude)[jall] / (denom * denom);
                dshift->element(j1i, j1i) -= (shift2) * lcovar * (*amplitude)[jall] / (denom * denom);
                dshift->element(j2i, j2i) += (shift2) * lcovar * (*amplitude)[jall] / (denom * denom);
                dshift->element(j3i, j3i) -= (shift2) * lcovar * (*amplitude)[jall] / (denom * denom);
              }
      }
    }
    {
      auto dtmp = den2_->copy();
      for (int i = 0; i != dtmp->ndim(); ++i)
        for (int j = 0; j != dtmp->mdim(); ++j)
          dtmp->element(j, i) += dshift->element(j, i);
      den2_ = dtmp;
    }
  }
#endif

  {
    // d_1^(2) -= <1|1><0|E_mn|0>     [Celani-Werner Eq. (A6)]
    auto dtmp = den2_->copy();
    for (int ist = 0; ist != nstates_; ++ist) {
      auto rdmtmp = rdm1all_->at(ist, ist)->matrix();
      const double factor = (*heff_)(ist, targetJ) * (*heff_)(ist, targetI);
      for (int i = nclosed; i != nclosed+nact; ++i)
        for (int j = nclosed; j != nclosed+nact; ++j) {
          dtmp->element(j, i) -= correlated_norm_lt_[ist] * (*rdmtmp)(j-nclosed, i-nclosed);
        }
    }
    dtmp->symmetrize();
    den2_ = dtmp;
  }

  if (info_->shift_imag()) {
    auto dtmp2 = den2_tt_->copy();
    for (int ist = 0; ist != nstates_; ++ist) {
      auto rdmtmp = rdm1all_->at(ist, ist)->matrix();
      const double factor = (*heff_)(ist, targetJ) * (*heff_)(ist, targetI);
      for (int i = nclosed; i != nclosed+nact; ++i)
        for (int j = nclosed; j != nclosed+nact; ++j) {
          dtmp2->element(j, i) -= correlated_norm_tt_[ist] * (*rdmtmp)(j-nclosed, i-nclosed) * factor;
        }
    }
    dtmp2->symmetrize();
    den2_tt_ = dtmp2;
  }

  shared_ptr<const Reference> ref = info_->ref();
  const MatView acoeff = coeff_->slice(nclosed+ncore, nclosed+ncore+nact);

  // code to calculate h+g(d). When add is false, h is not added
  auto focksub = [&](shared_ptr<const Matrix> moden, const MatView coeff, const bool add) {
    shared_ptr<const Matrix> jop = ref->geom()->df()->compute_Jop(make_shared<Matrix>(coeff * *moden ^ coeff));
    auto out = make_shared<Matrix>(acoeff % (add ? (*ref->hcore() + *jop) : *jop) * acoeff);
    shared_ptr<const DFFullDist> full = ref->geom()->df()->compute_half_transform(acoeff)->apply_J()->compute_second_transform(coeff)->swap();
    shared_ptr<DFFullDist> full2 = full->copy();
    full2 = full2->transform_occ1(moden);
    *out += *full->form_2index(full2, -0.5);
    return out;
  };

  if (!nocider) {
    shared_ptr<const Matrix> fock = focksub(ref->rdm1_mat(), coeff_->slice(0, ref->nocc()), true); // f
    shared_ptr<const Matrix> gd2 = focksub(den2_, coeff_->slice(ncore, coeff_->mdim()), false); // g(d2)

    // correct cideriv for fock derivative [Celani-Werner Eq. (C1), some terms in first and second lines]
    // y_I += (g[d^(2)]_ij - Nf_ij) <I|E_ij|0>
    for (int ist = 0; ist != nstates_; ++ist) {
      const Matrix op(*gd2 * (1.0/nstates_) - *fock * correlated_norm_lt_[ist]);
      shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          ci_deriv_->data(ist)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
    }

    if (info_->shift_imag()) {
      shared_ptr<const Matrix> gd2_tt = focksub(den2_tt_, coeff_->slice(ncore, coeff_->mdim()), false);
      for (int ist = 0; ist != nstates_; ++ist) {
        const double factor = (*heff_)(ist, targetJ) * (*heff_)(ist, targetI);
        const Matrix op2(*gd2_tt * (1.0/nstates_) - *fock * correlated_norm_tt_[ist] * factor);
        shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
        for (int i = 0; i != nact; ++i)
          for (int j = 0; j != nact; ++j)
            ci_deriv_->data(ist)->ax_plus_y(2.0*op2(j,i), deriv->data(j+i*nact));
      }
    }

    // y_I += <I|H|0> (for mixed states); taking advantage of the fact that unrotated CI vectors are eigenvectors
    if (targetJ == targetI) {
      // Gradient case. Special case of NACME with targetJ = targetI
      const Matrix ur(xmsmat_ ? *xmsmat_ * *heff_ : *heff_);
      const int target = targetJ;
      for (int ist = 0; ist != nstates_; ++ist)
        for (int jst = 0; jst != nstates_; ++jst)
          ci_deriv_->data(jst)->ax_plus_y(2.0*ur(ist,target)*(*heff_)(jst,target)*ref->energy(ist), info_orig_->ciwfn()->civectors()->data(ist));
    } else {
      // NACME case. targetJ and target I are separately used
      const Matrix ur(xmsmat_ ? *xmsmat_ * *heff_ : *heff_);
      for (int ist = 0; ist != nstates_; ++ist)
        for (int jst = 0; jst != nstates_; ++jst) {
          double urheff = (ur(ist,targetJ)*(*heff_)(jst,targetI) + ur(ist, targetI)*(*heff_)(jst,targetJ)) * ref->energy(ist);
          ci_deriv_->data(jst)->ax_plus_y(urheff, info_orig_->ciwfn()->civectors()->data(ist));
        }
    }

    // finally if this is XMS-CASPT2 gradient computation, we compute dcheck and contribution to y
    if (xmsmat_) {
      Matrix wmn(nstates_, nstates_);
      shared_ptr<Tensor> dc = rdm1_->clone();
      for (int i = 0; i != nstates_; ++i)
        for (int j = 0; j != i; ++j) {
          double cy = info_->ciwfn()->civectors()->data(j)->dot_product(ci_deriv_->data(i))
                    - info_->ciwfn()->civectors()->data(i)->dot_product(ci_deriv_->data(j));
          // If this is Full NACME, <U | dU/dX> contribution should be added
          if ((targetJ != targetI) && (nacmtype->is_full() || nacmtype->is_etf())) {
            cy += (pt2energy_[targetI] - pt2energy_[targetJ])
                * ((*heff_)(i,targetI) * (*heff_)(j,targetJ) - (*heff_)(j,targetI) * (*heff_)(i,targetJ));
          }
          wmn(j,i) = fabs(e0all_[j]-e0all_[i]) > 1.0e-12 ? -0.5 * cy / (e0all_[j]-e0all_[i]) : 0.0;
          wmn(i,j) = wmn(j,i);
          dc->ax_plus_y(wmn(j,i), rdm1all_->at(j, i));
          dc->ax_plus_y(wmn(i,j), rdm1all_->at(i, j));
        }
      dcheck_ = dc->matrix();

      // fill this into CI derivative. (Y contribution is done inside Z-CASSCF together with frozen core)
      shared_ptr<const Matrix> gdc = focksub(dcheck_, acoeff, false);
      for (int ist = 0; ist != nstates_; ++ist) {
        shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
        for (int jst = 0; jst != nstates_; ++jst) {
          Matrix op(*fock * wmn(jst, ist));
          if (ist == jst)
            op += *gdc * (1.0/nstates_) * 0.5;
          for (int i = 0; i != nact; ++i)
            for (int j = 0; j != nact; ++j)
              ci_deriv_->data(jst)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
        }
      }

      // also rotate cideriv back to the MS states
      btas::contract(1.0, *ci_deriv_->copy(), {0,1,2}, (*xmsmat_), {3,2}, 0.0, *ci_deriv_, {0,1,3});
    }
  }

  if (info_->shift_imag()) {
    auto dtmp = den2_->copy();
    dtmp->ax_plus_y(1.0, den2_tt_);
    den2_ = dtmp;
  }

  // restore original energy
  energy_ = pt2energy_;
  timer.tick_print("Postprocessing SMITH");
}


#endif
