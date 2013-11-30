//
// BAGEL - Parallel electron correlation program.
// Filename: ras/dist_form_sigma.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <map>
#include <algorithm>

#include <src/ras/dist_form_sigma.h>
#include <src/math/sparsematrix.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Olsen */
// This is how this is accessed from RASCI
shared_ptr<DistRASDvec> DistFormSigmaRAS::operator()(shared_ptr<const DistRASDvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();

  const int norb = det->norb();
  auto g = make_shared<Matrix>(norb, norb);
  for (int k = 0, kl = 0; k < norb; ++k) {
    for (int l = 0; l < k; ++l, ++kl) {
      { // g_kl
        double val = jop->mo1e(kl) - jop->mo2e_hz(k, k, k, l);
        for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,l);
        g->element(l, k) = val;
      }

      { // g_lk
        double val = jop->mo1e(kl);
        for (int j = 0; j < l; ++j) val -= jop->mo2e_hz(l,j,j,k);
        g->element(k, l) = val;
      }
    }
    // g_kk
    double val = jop->mo1e(kl) - 0.5*jop->mo2e_hz(k,k,k,k);
    for (int j = 0; j < k; ++j) val -= jop->mo2e_hz(k,j,j,k);
    g->element(k, k) = val;
    ++kl;
  }

  auto sigmavec = make_shared<DistRASDvec>(det, nstate);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(0);
    if (conv[istate]) continue;
    shared_ptr<const DistRASCivec> cc = ccvec->data(istate);
    shared_ptr<DistRASCivec> sigma = sigmavec->data(istate);

    // start transpose
    shared_ptr<DistRASCivec> cctrans = cc->transpose();
    pdebug.tick_print("transpose");

    // (taskbb)
    sigma_bb(cc, sigma, g->data(), jop->mo2e_ptr());
    pdebug.tick_print("taskbb");

    // finish transpose
    cctrans->transpose_wait();
    pdebug.tick_print("wait");

    // (taskaa)
    shared_ptr<DistRASCivec> strans = cctrans->clone();
    sigma_bb(cctrans, strans, g->data(), jop->mo2e_ptr());
    pdebug.tick_print("taskaa");

    // start transpose back
    shared_ptr<DistRASCivec> saa = strans->transpose(det);
    pdebug.tick_print("transpose1");

    // (taskab) alpha-beta contributions
    sigma_ab(cc, sigma, jop->mo2e_ptr());
    pdebug.tick_print("taskab");

    // finish transpose back
    saa->transpose_wait();
    pdebug.tick_print("wait1");
    sigma->ax_plus_y(1.0, *saa);
  }
  return sigmavec;
}

// A bit of a temporary hack for 1e terms
shared_ptr<DistRASDvec> DistFormSigmaRAS::operator()(shared_ptr<const DistRASDvec> ccvec, const double* mo1e) const {
  const int nstate = ccvec->ij();
  shared_ptr<const RASDeterminants> det = ccvec->det();
  const int norb = det->norb();

  auto sigmavec = make_shared<DistRASDvec>(det, nstate);

  unique_ptr<double[]> blank2e(new double[norb*norb*norb*norb]);
  fill_n(blank2e.get(), norb*norb*norb*norb, 0.0);

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<const DistRASCivec> cc = ccvec->data(istate);
    shared_ptr<DistRASCivec> sigma = sigmavec->data(istate);

    // start transpose
    shared_ptr<DistRASCivec> cctrans = cc->transpose();

    // (taskbb)
    sigma_bb(cc, sigma, mo1e, blank2e.get());

    shared_ptr<DistRASCivec> strans = cctrans->clone();
    cctrans->transpose_wait();

    // (taskaa)
    sigma_bb(cctrans, strans, mo1e, blank2e.get());

    shared_ptr<DistRASCivec> saa = strans->transpose();
    saa->transpose_wait();

    sigma->ax_plus_y(1.0, *saa);
  }

  return sigmavec;
}

// sigma_bb requires no communication
void DistFormSigmaRAS::sigma_bb(shared_ptr<const DistRASCivec> cc, shared_ptr<DistRASCivec> sigma, const double* g, const double* mo2e) const {
  shared_ptr<const RASDeterminants> det = cc->det();
  assert(*det == *sigma->det());

  const int norb = det->norb();
  const size_t lb = det->lenb();

  for (auto& spaceiter : det->stringspaceb()) {
    shared_ptr<const StringSpace> ispace = spaceiter.second;
    auto F = make_shared<Matrix>(lb, ispace->size());
    double* fdata = F->data();
    for (size_t ib = 0; ib < ispace->size(); ++ib, fdata+=lb) {
      for (auto& iterkl : det->phib(ib + ispace->offset())) {
        fdata[iterkl.source] += static_cast<double>(iterkl.sign) * g[iterkl.ij];
        for (auto& iterij : det->phib(iterkl.source)) {
          if (iterij.ij < iterkl.ij) continue;
          const int ii = iterij.ij/norb;
          const int jj = iterij.ij%norb;
          const int kk = iterkl.ij/norb;
          const int ll = iterkl.ij%norb;
          fdata[iterij.source] += static_cast<double>(iterkl.sign * iterij.sign) * (iterkl.ij == iterij.ij ? 0.5 : 1.0) * mo2e[ii + kk*norb + norb*norb*(jj + ll * norb)];
        }
      }
    }

    // F is finished, apply to the right places
    for (auto& iblock : cc->blocks()) {
      if (!iblock) continue;
      if (!det->allowed(ispace, iblock->stringa())) continue;
      shared_ptr<DistRASBlock<double>> tblock = sigma->block(ispace, iblock->stringa());

      dgemm_("T", "N", ispace->size(), tblock->asize(), iblock->lenb(), 1.0,
        F->element_ptr(iblock->stringb()->offset(), 0), F->ndim(), iblock->local(), iblock->lenb(), 1.0, tblock->local(), tblock->lenb());
    }
  }
}

void DistFormSigmaRAS::sigma_ab(shared_ptr<const DistRASCivec> cc, shared_ptr<DistRASCivec> sigma, const double* mo2e) const {
  assert(cc->det() == sigma->det());
  shared_ptr<const RASDeterminants> det = cc->det();

  const int norb = det->norb();
  list<tuple<int, shared_ptr<Matrix>, shared_ptr<Matrix>, shared_ptr<const StringSpace>, int>> requests;

  // mapping space offsets to process bounds
  map<size_t, tuple<size_t, size_t>> bounds_map;
  map<size_t, vector<int>> scattering_map;
  for (auto& spaceiter : det->stringspacea()) {
    shared_ptr<const StringSpace> sp = spaceiter.second;
    StaticDist d(sp->size(), mpi__->size());
    bounds_map.emplace(sp->offset(), d.range(mpi__->rank()));
    vector<int> scat(mpi__->size());
    for (int i = 0; i < mpi__->size(); ++i)
      scat[i] = d.size(i);
    scattering_map.emplace(sp->offset(), scat);
  }


  map<size_t, map<size_t, pair<vector<tuple<size_t, int, size_t>>, shared_ptr<SparseMatrix>>>> Fmatrices;

  // Builds prototypes for the sparse F matrices as well as vectors that contain all of the information necessary to update the matrices as they are needed
  const int nspaces = det->stringspacea().size();
  for (auto& spaceiter : det->stringspacea()) {
    shared_ptr<const StringSpace> ispace = spaceiter.second;
    const size_t la = ispace->size();

    // These are for building the initial versions of the sparse matrices
    vector<vector<double>> data(nspaces);
    vector<vector<int>> cols(nspaces);
    vector<vector<int>> rind(nspaces);
    vector<vector<tuple<size_t, int, size_t>>> sparse_info(nspaces);

    vector<pair<size_t, size_t>> bounds;
    for (auto& isp : bounds_map)
      bounds.emplace_back(isp.first + get<0>(isp.second), isp.first + get<1>(isp.second));
    assert(bounds.size() == nspaces);

    for (int ia = 0; ia < la; ++ia) {
      map<size_t, double> row;
      map<size_t, vector<tuple<int, int>>> row_positions;
      for (auto& iter : det->phia(ia + ispace->offset())) {
        const int kk = iter.ij/norb;
        const int ll = iter.ij%norb;
        row_positions[iter.source].emplace_back(iter.sign, ll*norb + kk*norb*norb*norb);
      }

      for (int sp = 0; sp < nspaces; ++sp) rind[sp].push_back(data[sp].size() + 1);

      auto ibound = bounds.begin();
      for (auto& irow : row_positions) {
        for (int isp = 0; isp < nspaces; ++isp) {
          if ( irow.first >= bounds[isp].first && irow.first < bounds[isp].second) {
            const size_t pos = data[isp].size();
            for (auto& i : irow.second)
              sparse_info[isp].emplace_back(pos, get<0>(i), get<1>(i));
            cols[isp].push_back(irow.first + 1 - bounds[isp].first);
            data[isp].push_back(1.0);
            continue;
          }
        }
      }
    }

    // Map connects the offset of a stringspace to a corresponding SparseMatrix and all the info needed to reconstruct it
    map<size_t, pair<vector<tuple<size_t, int, size_t>>, shared_ptr<SparseMatrix>>> Fmap;

    // TODO is using the starting bound okay here?
    for (int isp = 0; isp < nspaces; ++isp) {
      if (data[isp].size() > 0) {
        rind[isp].push_back(data[isp].size() + 1);
        const int mdim = bounds.at(isp).second - bounds.at(isp).first;
        Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), make_shared<SparseMatrix>(la, mdim, data[isp], cols[isp], rind[isp])));
      }
      else Fmap.emplace(bounds[isp].first, make_pair(move(sparse_info[isp]), shared_ptr<SparseMatrix>()));
    }

    Fmatrices.emplace(ispace->offset(), move(Fmap));
  }

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      // L(I), R(I), sign(I) building
      const size_t phisize = accumulate(det->phib_ij(ij).begin(), det->phib_ij(ij).end(), 0ull, [] (size_t i, const RAS::DMapBlock& m) { return i + m.size(); });
      if (phisize == 0) continue;

      map<pair<size_t, size_t>, shared_ptr<Matrix>> Cp_map;

      // gathering
      for ( auto& iphiblock : det->phib_ij(ij) ) {
        vector<shared_ptr<const DistRASBlock<double>>> blks = cc->allowed_blocks<1>(iphiblock.space());
        for (auto& iblock : blks) {
          auto tmp = make_shared<Matrix>(iblock->asize(), iphiblock.size());
          double* targetdata = tmp->data();

          const size_t lb = iblock->lenb();

          for (auto& iphi : iphiblock) {
            double sign = static_cast<double>(iphi.sign);
            const double* sourcedata = iblock->local() + iphi.source;

            for (size_t i = 0; i < iblock->asize(); ++i, ++targetdata, sourcedata+=lb)
              *targetdata = *sourcedata * sign;
          }

          // TODO is using the offset okay here?
          Cp_map.emplace(make_pair(iblock->stringa()->offset(), iphiblock.offset()), tmp);
        }
      }

      // build V(I), block by block
      for (auto& spaceiter : det->stringspacea()) {
        shared_ptr<const StringSpace> ispace = spaceiter.second;
        const size_t la = ispace->size();

        // Beware, this COULD be a memory problem
        auto Vt = make_shared<Matrix>(la, phisize);
        auto& Fmap = Fmatrices.at(ispace->offset());

        // Replace data in SparseMatrix
        const double* mo2e_ij = mo2e + i + norb*norb*j;
        for (auto& f : Fmap) {
          shared_ptr<SparseMatrix> sparse = f.second.second;
          if (sparse) {
            sparse->zero();
            double* fdata = sparse->data();
            vector<tuple<size_t, int, size_t>>& replace_data = f.second.first;
            for_each(replace_data.begin(), replace_data.end(), [&fdata, &mo2e_ij] (const tuple<size_t, int, size_t>& i)
              { fdata[get<0>(i)] += static_cast<double>(get<1>(i)) * mo2e_ij[get<2>(i)]; });
          }
        }

        for (auto& iphiblock : det->phib_ij(ij)) {
          vector<shared_ptr<const StringSpace>> allowed_spaces = det->allowed_spaces<1>(iphiblock.space());
          for (auto& mult_space : allowed_spaces) {
            shared_ptr<Matrix> Cp_block = Cp_map.at(make_pair(mult_space->offset(), iphiblock.offset()));
            shared_ptr<SparseMatrix> Ft_block = Fmap[mult_space->offset() + get<0>(bounds_map[mult_space->offset()])].second;
            if (Ft_block) {
              auto Vt_block = make_shared<Matrix>(*Ft_block * *Cp_block);
              Vt->add_block(1.0, 0, iphiblock.offset(), Vt_block->ndim(), Vt_block->mdim(), Vt_block);
            }
          }
        }

        // Add up contributions from each node
        const size_t astart = get<0>(bounds_map[ispace->offset()]);
        const size_t aend = get<1>(bounds_map[ispace->offset()]);
        const size_t asize = aend - astart;

        shared_ptr<Matrix> V = Vt->transpose();

        vector<int> recv_counts(mpi__->size());
        vector<int>& asizes = scattering_map[ispace->offset()];
        for (int i = 0; i < mpi__->size(); ++i)
          recv_counts[i] = asizes[i] * phisize;

        auto V_chunk = make_shared<Matrix>(phisize, aend - astart);
        //mpi__->reduce_scatter(V->data(), V_chunk->data(), recv_counts.data());
        const int rq = mpi__->ireduce_scatter(V->data(), V_chunk->data(), recv_counts.data());
        requests.emplace_back(rq, V, V_chunk, ispace, ij);
      }

        // scatter
      for (auto ir = requests.begin(); ir != requests.end(); ) {
        const int rq = get<0>(*ir);
        mpi__->wait(rq);
        {
        //if (mpi__->test(rq)) {
          shared_ptr<Matrix> Vt_chunk = get<2>(*ir)->transpose();
          shared_ptr<const StringSpace> ispace = get<3>(*ir);
          const int ij = get<4>(*ir);

          const size_t astart = get<0>(bounds_map[ispace->offset()]);
          const size_t aend = get<1>(bounds_map[ispace->offset()]);
          const size_t asize = aend - astart;
          if (asize == 0) continue;

          const double* vdata = Vt_chunk->data();
          for (auto& iphiblock : det->phib_ij(ij) ) {
            for (auto& iphi : iphiblock) {
              shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(iphi.target));
              if (det->allowed(ispace, betaspace)) {

                shared_ptr<DistRASBlock<double>> sgblock = sigma->block(betaspace, ispace);
                double* const targetdata = sgblock->local() + iphi.target - betaspace->offset();

                const size_t lb = sgblock->lenb();

                for (size_t ii = sgblock->astart(), jj = 0; ii < sgblock->aend(); ++ii, ++jj)
                  targetdata[jj*lb] += vdata[jj];
              }

              vdata += asize;
            }
          }
          ir = requests.erase(ir);
        }
        //else {
        //  ++ir;
        //}
      }
    }
  }

  bool done;
  do {
    done = true;
    for (auto ir = requests.begin(); ir != requests.end(); ) {
      const int rq = get<0>(*ir);
      mpi__->wait(rq);
      {
      //if (mpi__->test(rq)) {
        shared_ptr<Matrix> Vt_chunk = get<2>(*ir)->transpose();
        shared_ptr<const StringSpace> ispace = get<3>(*ir);
        const int ij = get<4>(*ir);

        const size_t astart = get<0>(bounds_map[ispace->offset()]);
        const size_t aend = get<1>(bounds_map[ispace->offset()]);
        const size_t asize = aend - astart;

        const double* vdata = Vt_chunk->data();
        for (auto& iphiblock : det->phib_ij(ij) ) {
          for (auto& iphi : iphiblock) {
            shared_ptr<const StringSpace> betaspace = det->space<1>(det->stringb(iphi.target));
            if (det->allowed(ispace, betaspace)) {

              shared_ptr<DistRASBlock<double>> sgblock = sigma->block(betaspace, ispace);
              double* const targetdata = sgblock->local() + iphi.target - betaspace->offset();

              const size_t lb = sgblock->lenb();

              for (size_t ii = sgblock->astart(), jj = 0; ii < sgblock->aend(); ++ii, ++jj)
                targetdata[jj*lb] += vdata[jj];
            }

            vdata += asize;
          }
        }
        ir = requests.erase(ir);
      }
      //else {
      //  ++ir;
      //  done = false;
      //}
    }
    if (!done) std::this_thread::sleep_for(sleeptime__);
  } while (!done);
  mpi__->barrier();
}
