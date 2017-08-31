//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diracgrad.cc
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

#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/alpha.h>
#include <src/df/reldffull.h>
#include <src/scf/dhf/dfock.h>
#include <src/wfn/relreference.h>
#include <src/integral/rys/gsmallnaibatch.h>

using namespace std;
using namespace bagel;

#define LOCAL_TIMING

template<>
shared_ptr<GradFile> GradEval<Dirac>::compute(const string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  assert(ref_->nact() == 0);
  geom_ = task_->geom();

  Timer timer;
#ifdef LOCAL_TIMING
  Timer ptime(0);
#endif
  // density matrix
  shared_ptr<const RelReference> ref = dynamic_pointer_cast<const RelReference>(ref_);
  shared_ptr<const ZMatrix> coeff = ref->relcoeff()->slice_copy(0, 2*ref->nocc());
  auto den = make_shared<const ZMatrix>(*coeff ^ *coeff);
  if (ref->gaunt())
    throw runtime_error("Dirac--Hartree--Fock gradients are currently unavailable with the Gaunt or full Breit interaction.");

  // energy-weighted density matrix
  shared_ptr<ZMatrix> ecoeff = coeff->copy();
  const VectorB& eig = ref->eig();
  for (int i = 0; i != 2*ref->nocc(); ++i)
    zscal_(ecoeff->ndim(), eig(i), ecoeff->element_ptr(0, i), 1);
  auto eden = make_shared<const ZMatrix>(*coeff ^ *ecoeff);

  const int nbasis = geom_->nbasis();

  // NAI density (L+L+, L-L-)
  shared_ptr<ZMatrix> nden  =  den->get_submatrix(0, 0, nbasis, nbasis);
                     *nden += *den->get_submatrix(nbasis, nbasis, nbasis, nbasis);
  // kinetic density [den] 2*(S+L+, S-L-) - (S+S+, S-S-); [eden] -(S+S+, S-S-)/2c^2
  shared_ptr<ZMatrix> kden  =  den->get_submatrix(0, 2*nbasis, nbasis, nbasis);
                     *kden += *den->get_submatrix(nbasis, 3*nbasis, nbasis, nbasis);
                     *kden *= complex<double>(2.0);
                     *kden -= *den->get_submatrix(2*nbasis, 2*nbasis, nbasis, nbasis);
                     *kden -= *den->get_submatrix(3*nbasis, 3*nbasis, nbasis, nbasis);
  shared_ptr<ZMatrix> lden  = eden->get_submatrix(2*nbasis, 2*nbasis, nbasis, nbasis);
                     *lden +=*eden->get_submatrix(3*nbasis, 3*nbasis, nbasis, nbasis);
                     *lden /= complex<double>(2.0*pow(c__,2));
                     *kden -= *lden;
  // overlap density
  shared_ptr<ZMatrix> sden  = eden->get_submatrix(0, 0, nbasis, nbasis);
                     *sden +=*eden->get_submatrix(nbasis, nbasis, nbasis, nbasis);

  // nden, kden, sden (minus sign is taken care of inside)
  vector<shared_ptr<GradTask>> task = contract_grad1e<GradTask1>(nden->get_real_part(), kden->get_real_part(), sden->get_real_part());

  // small NAI part..
  map<int, shared_ptr<Sigma>> sigma;
  sigma.emplace(Comp::X, make_shared<Sigma>(Comp::X));
  sigma.emplace(Comp::Y, make_shared<Sigma>(Comp::Y));
  sigma.emplace(Comp::Z, make_shared<Sigma>(Comp::Z));
  auto sp = make_shared<ZMatrix>(4,1,true); sp->element(2,0) = 1;
  auto sm = make_shared<ZMatrix>(4,1,true); sm->element(3,0) = 1;
  map<int, shared_ptr<ZMatrix>> XY{ make_pair(2, sp), make_pair(3, sm) };

  // target data area
  vector<int> xyz{Comp::X, Comp::Y, Comp::Z};
  map<pair<int,int>, shared_ptr<ZMatrix>> mat;
  for (auto& i : xyz)
    for (auto& j : xyz)
      if (i <= j)
        mat.emplace(make_pair(i,j), make_shared<ZMatrix>(nbasis, nbasis));

  for (auto& s0 : XY) { // bra
    for (auto& s1 : XY) { // ket
      shared_ptr<ZMatrix> data = den->get_submatrix(s0.first*nbasis, s1.first*nbasis, nbasis, nbasis);
      for (auto& w0 : sigma) {
        for (auto& w1 : sigma) {
          const complex<double> c = ((*w0.second * *s0.second) % (*w1.second * *s1.second)).element(0,0);
          const int small = min(w0.first, w1.first);
          const int large = max(w0.first, w1.first);
          mat[{small, large}]->ax_plus_y(c, data);
        }
      }
    }
  }
  assert(mat.size() == 6);
  array<shared_ptr<const Matrix>,6> rmat;
  auto riter = rmat.begin();
  for (auto& i : mat) *riter++ = i.second->get_real_part();

  // *** adding task here ****
  {
    vector<shared_ptr<GradTask>> tmp = contract_gradsmall1e(rmat);
    task.insert(task.end(), tmp.begin(), tmp.end());
  }

  if (geom_->has_finite_nucleus()) {
    vector<shared_ptr<GradTask>> tmp = contract_grad1e_fnai(rmat);
    task.insert(task.end(), tmp.begin(), tmp.end());
  }
#ifdef LOCAL_TIMING
  mpi__->barrier();
  ptime.tick_print("Onebody part");
#endif


  // two-electron contributions.
  {
    // make blocks of coefficients
    array<shared_ptr<const Matrix>, 4> rocoeff;
    array<shared_ptr<const Matrix>, 4> iocoeff;
    array<shared_ptr<const Matrix>, 4> trocoeff;
    array<shared_ptr<const Matrix>, 4> tiocoeff;
    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> ocoeff = coeff->get_submatrix(i*geom_->nbasis(), 0, geom_->nbasis(), 2*ref->nocc());
      rocoeff[i] = ocoeff->get_real_part();
      iocoeff[i] = ocoeff->get_imag_part();
      trocoeff[i] = rocoeff[i]->transpose();
      tiocoeff[i] = iocoeff[i]->transpose();
    }
    // (0) get AO integras
    // get individual df dist objects for each block and add df to dfs
    list<shared_ptr<RelDFHalf>> half_complex_exch;

    const bool external_half = !task_->half().empty();
    if (external_half) {
      // (1-3) Reuse half-transform integrals if possible
      half_complex_exch = task_->half();
      task_->discard_half();
    } else {
      vector<shared_ptr<const DFDist>> dfs = geom_->dfs()->split_blocks();
      dfs.push_back(geom_->df());

      // (1) make RelDF objects from AO integrals
      list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

      // (2) first-transform
      list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, coeff->slice_copy(0, 2*ref->nocc()));

      // (3) split and factorize
      for (auto& i : half_complex) {
        list<shared_ptr<RelDFHalf>> tmp = i->split(false);
        half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
      }
      half_complex.clear();
      DFock::factorize(half_complex_exch);
    }
#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("first transformed");
#endif

    // (4) compute C matrix
    shared_ptr<RelCDMatrix> cd;
    for (auto& j : half_complex_exch) {
      for (auto& i : j->basis()) {
        if (cd) {
          *cd += RelCDMatrix(j, i, trocoeff, tiocoeff, geom_->df()->data2(), external_half ? 1 : 2);
        } else {
          cd = make_shared<RelCDMatrix>(j, i, trocoeff, tiocoeff, geom_->df()->data2(), external_half ? 1 : 2);
        }
      }
    }

    // (5) compute (gamma|ij)
    list<shared_ptr<RelDFFull>> dffull;
    for (auto& i : half_complex_exch) {
      auto tmp = make_shared<RelDFFull>(i, rocoeff, iocoeff);
      dffull.push_back(external_half ? tmp->apply_J() : tmp->apply_JJ());
    }
    DFock::factorize(dffull);
    dffull.front()->scale(dffull.front()->fac()); // take care of the factor
    assert(dffull.size() == 1);
#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("second transformed");
#endif

    // (6) two-index gamma
    shared_ptr<VectorB> cdr = cd->get_real_part();
    assert(cd->get_imag_part()->rms() < 1.0e-10); // by symmetry the imaginary part is zero
    auto gamma2 = make_shared<const Matrix>((*cdr ^ *cdr) - *dffull.front()->form_aux_2index_real());

#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("gamma2");
#endif

    // *** adding task here ****
    {
      vector<shared_ptr<GradTask>> task2 = contract_grad2e_2index(gamma2);
      task.insert(task.end(), task2.begin(), task2.end());
    }

    // (7) first back transformation (gamma|is^Y)
    list<shared_ptr<RelDFHalfB>> dfhalfb = dffull.front()->back_transform(rocoeff, iocoeff);

#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("first backtransformed");
#endif

    // (8) second back transformation (gamma|r^Xs^Y) and immediately rearrange to (gamma|r^w s^Y)
    map<pair<int,int>,shared_ptr<DFDist>> gamma3;
    for (auto& half : dfhalfb) {
      const int cbasis = half->basis();
      if (cbasis == Basis::LP || cbasis == Basis::LM) {
        // large component
        pair<int,int> key {Comp::L, Comp::L};
        auto iter = gamma3.find(key);
        if (iter == gamma3.end()) {
          gamma3.emplace(key, half->back_transform(rocoeff[cbasis], iocoeff[cbasis]));
        } else {
          iter->second->ax_plus_y(1.0, half->back_transform(rocoeff[cbasis], iocoeff[cbasis])); // TODO redundant copy, but probably fine
        }
      } else {
        array<int,2> SS {{ Basis::SP, Basis::SM }};
        for (auto& cbasis1 : SS) {
          shared_ptr<DFDist> rdf = half->back_transform(rocoeff[cbasis1], iocoeff[cbasis1]);
          shared_ptr<DFDist> idf = half->back_transform(rocoeff[cbasis1], iocoeff[cbasis1], true);
          assert(XY.find(cbasis) != XY.end() && XY.find(cbasis1) != XY.end());
          auto s0 = XY[cbasis];
          auto s1 = XY[cbasis1];
          for (auto& w0 : sigma) {
            for (auto& w1 : sigma) {
              if (w0.first <= w1.first) {
                // calculate k^ww'_XY
                const complex<double> kwwxx = conj(((*w0.second * *s0) % (*w1.second * *s1)).element(0,0));
                const bool imag = fabs(kwwxx.imag()) > 1.0e-20;
                assert(!imag || fabs(kwwxx.real()) < 1.0e-20);
                const double fac = 1.0 * (imag ? -kwwxx.imag() : kwwxx.real()); // -1 comes from the prefactor of exchange

                // TODO we can use symmetry
                const int small = min(w0.first, w1.first);
                const int large = max(w0.first, w1.first);

                pair<int,int> key {small, large};
                auto iter = gamma3.find(key);
                if (iter == gamma3.end()) {
                  gamma3.emplace(key, imag ? idf->copy() : rdf->copy());
                  gamma3[key]->scale(fac);
                } else {
                  iter->second->ax_plus_y(fac, imag ? idf : rdf);
                }
              }
            }
          }
        }
      }
    }

    // (9) direct product contributions
    map<pair<int,int>,shared_ptr<const Matrix>> wden;
    for (auto& r : mat)
      wden.emplace(r.first, r.second->get_real_part());
    wden.emplace(make_pair(Comp::L,Comp::L), nden->get_real_part()); // large-large case
    for (auto& w : wden) {
      auto iter = gamma3.find(w.first);
      assert(iter != gamma3.end());
      iter->second->add_direct_product(cdr, w.second, -1.0);
    }

    // minus one
    for (auto& i : gamma3) {
      i.second->scale(-1.0);
    }

#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("second backtransformed");
#endif

    // *** adding task here ****
    { // large-large
      auto iter = gamma3.find({Comp::L,Comp::L});
      assert(iter != gamma3.end());
      // transform to the shell-boundary format
      iter->second->shell_boundary_3index();
      vector<shared_ptr<GradTask>> task3 = contract_grad2e(iter->second);
      task.insert(task.end(), task3.begin(), task3.end());
    }
    { // small-small
      array<shared_ptr<const DFDist>,6> gs;
      int icnt = 0;
      for (auto& i : xyz)
        for (auto& j : xyz)
          if (i <= j) {
            auto iter = gamma3.find({i,j});
            assert(iter != gamma3.end());
            // transform to the shell-boundary format
            iter->second->shell_boundary_3index();
            gs[icnt++] = iter->second;
          }
      vector<shared_ptr<GradTask>> task3 = contract_grad2e(gs);
      task.insert(task.end(), task3.begin(), task3.end());
    }
  }

  // compute
  TaskQueue<shared_ptr<GradTask>> tq(move(task));
  tq.compute();

  // adds nuclear contributions
  *grad_ += *geom_->compute_grad_vnuc();

  // allreduce
  grad_->allreduce();

#ifdef LOCAL_TIMING
    mpi__->barrier();
    ptime.tick_print("gradient computation");
#endif

  grad_->print();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  energy_ = ref_->energy(0);

  return grad_;
}

