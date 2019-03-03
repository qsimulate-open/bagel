//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dftgrid.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#include <numeric>
#include <src/scf/ks/dftgrid.h>
#include <src/scf/ks/lebedevlist.h>
#include <src/scf/ks/xcfunc.h>
#include <src/util/f77.h>
#include <src/util/constants.h>
#include <src/util/taskqueue.h>
#include <src/util/timer.h>
#include <src/util/parallel/staticdist.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

const static LebedevList lebedev;


vector<shared_ptr<const Matrix>> DFTGrid_base::compute_rho_sigma(shared_ptr<const XCFunc> func, shared_ptr<const Matrix> mat,
                                                         unique_ptr<double[]>& rho, unique_ptr<double[]>& sigma,
                                                         unique_ptr<double[]>& rhox, unique_ptr<double[]>& rhoy, unique_ptr<double[]>& rhoz) const {
  vector<shared_ptr<const Matrix>> out;
  auto orb = make_shared<Matrix>(*mat % *grid_->basis());
  if (func->lda()) {
    assert(orb->mdim() == grid_->size());
    for (size_t i = 0; i != orb->mdim(); ++i) {
      rho[i] = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orb->element_ptr(0, i), 1);
    }
    out = vector<shared_ptr<const Matrix>>{orb};
  } else {
    auto orbx = make_shared<Matrix>(*mat % *grid_->gradx());
    auto orby = make_shared<Matrix>(*mat % *grid_->grady());
    auto orbz = make_shared<Matrix>(*mat % *grid_->gradz());
    for (size_t i = 0; i != orb->mdim(); ++i) {
      rho[i] = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orb->element_ptr(0, i), 1);
      const double sigx = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orbx->element_ptr(0, i), 1);
      const double sigy = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orby->element_ptr(0, i), 1);
      const double sigz = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orbz->element_ptr(0, i), 1);
      sigma[i] = 4*(sigx*sigx + sigy*sigy + sigz*sigz);
      rhox[i] = 2*sigx;
      rhoy[i] = 2*sigy;
      rhoz[i] = 2*sigz;
    }
    out = vector<shared_ptr<const Matrix>>{orb, orbx, orby, orbz};
  }
  return out;
}


namespace bagel {
class ExcVxcTask {
  protected:
    const size_t size;
    const double* rho;
    const double* sigma;
    double* exc;
    double* vxc;
    double* vxc2;
    shared_ptr<const XCFunc> func;
  public:
    ExcVxcTask(const size_t n, const double* r, const double* s, double* e, double* v, double* v2, shared_ptr<const XCFunc> f)
     : size(n), rho(r), sigma(s), exc(e), vxc(v), vxc2(v2), func(f) { }
    void compute() {
      func->compute_exc_vxc(size, rho, sigma, exc, vxc, vxc2);
    }
};
class VxcTask {
  protected:
    const size_t size;
    const double* rho;
    const double* sigma;
    double* vxc;
    double* vxc2;
    shared_ptr<const XCFunc> func;
  public:
    VxcTask(const size_t n, const double* r, const double* s, double* v, double* v2, shared_ptr<const XCFunc> f)
     : size(n), rho(r), sigma(s), vxc(v), vxc2(v2), func(f) { }
    void compute() {
      func->compute_vxc(size, rho, sigma, vxc, vxc2);
    }
};
}


tuple<shared_ptr<const Matrix>,double> DFTGrid_base::compute_xc(shared_ptr<const XCFunc> func, shared_ptr<const Matrix> mat) const {
  Timer time;

  unique_ptr<double[]> rho(new double[grid_->size()]);
  unique_ptr<double[]> sigma, rhox, rhoy, rhoz;
  if (!func->lda()) {
    sigma = unique_ptr<double[]>(new double[grid_->size()]);
    rhox  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoy  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoz  = unique_ptr<double[]>(new double[grid_->size()]);
  }

  compute_rho_sigma(func, mat, rho, sigma, rhox, rhoy, rhoz);

  time.tick_print("rho, sigma");

  unique_ptr<double[]> exc(new double[grid_->size()]);
  unique_ptr<double[]> vxc(new double[grid_->size()*(func->lda()?1:2)]);

  StaticDist dist(grid_->size(), min(resources__->max_num_threads()*100, grid_->size()));
  vector<pair<size_t, size_t>> table = dist.atable();

  TaskQueue<ExcVxcTask> tasks(table.size());
  for (auto& i : table) {
    const size_t n = i.first;
    tasks.emplace_back(i.second, rho.get()+n, (!func->lda() ? sigma.get()+n : nullptr),
                               exc.get()+n, vxc.get()+n, (!func->lda() ? vxc.get()+n+grid_->size() : nullptr), func);
  }
  tasks.compute();
  time.tick_print("exc+vxc");

  auto out = make_shared<Matrix>(mol_->nbasis(), mol_->nbasis());
  double en = 0.0;

  auto scal = make_shared<Matrix>(mol_->nbasis(), grid_->size());
  if (func->lda()) {
    for (size_t i = 0; i != scal->mdim(); ++i) {
      daxpy_(scal->ndim(), vxc[i]*grid_->weight(i), grid_->basis()->element_ptr(0, i), 1, scal->element_ptr(0, i), 1);
      en += exc[i] * rho[i] * grid_->weight(i);
    }
  } else {
    for (size_t i = 0; i != scal->mdim(); ++i) {
      daxpy_(scal->ndim(), vxc[i]*grid_->weight(i), grid_->basis()->element_ptr(0, i), 1, scal->element_ptr(0, i), 1);
      daxpy_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhox[i], grid_->gradx()->element_ptr(0, i), 1, scal->element_ptr(0, i), 1);
      daxpy_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhoy[i], grid_->grady()->element_ptr(0, i), 1, scal->element_ptr(0, i), 1);
      daxpy_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhoz[i], grid_->gradz()->element_ptr(0, i), 1, scal->element_ptr(0, i), 1);
      en += exc[i] * rho[i] * grid_->weight(i);
    }
  }
  *out += *scal ^ *grid_->basis();
  out->symmetrize();

  time.tick_print("contraction");
  return make_tuple(out, en);
}


shared_ptr<const GradFile> DFTGrid_base::compute_xcgrad(shared_ptr<const XCFunc> func, shared_ptr<const Matrix> mat) const {
  auto out = make_shared<GradFile>(mol_->natom());

  unique_ptr<double[]> rho(new double[grid_->size()]);
  unique_ptr<double[]> sigma, rhox, rhoy, rhoz;
  if (!func->lda()) {
    sigma = unique_ptr<double[]>(new double[grid_->size()]);
    rhox  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoy  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoz  = unique_ptr<double[]>(new double[grid_->size()]);
  }
  unique_ptr<double[]> vxc(new double[grid_->size()*(func->lda()?1:2)]);

  vector<shared_ptr<const Matrix>> orb = compute_rho_sigma(func, mat, rho, sigma, rhox, rhoy, rhoz);

  StaticDist dist(grid_->size(), min(resources__->max_num_threads()*100, grid_->size()));
  vector<pair<size_t, size_t>> table = dist.atable();

  TaskQueue<VxcTask> tasks(table.size());
  for (auto& i : table) {
    const size_t n = i.first;
    tasks.emplace_back(i.second, rho.get()+n, (!func->lda() ? sigma.get()+n : nullptr),
                            vxc.get()+n, (!func->lda() ? vxc.get()+n+grid_->size() : nullptr), func);
  }
  tasks.compute();

  // in GGA, we need nabla^2 basis // TODO I guess this should be more efficient..
  array<shared_ptr<Matrix>,6> grad2;
  if (!func->lda())
    grad2 = grid_->compute_grad2();

  // loop over target atom
  size_t offset = 0;
  int n = 0;
  for (auto& b : mol_->atoms()) {
    shared_ptr<const Matrix> bmat = mat->cut(offset, offset+b->nbasis());
    array<shared_ptr<const Matrix>,3> d1mat;
    d1mat[0] = make_shared<const Matrix>(*bmat % *grid_->gradx()->cut(offset, offset+b->nbasis()));
    d1mat[1] = make_shared<const Matrix>(*bmat % *grid_->grady()->cut(offset, offset+b->nbasis()));
    d1mat[2] = make_shared<const Matrix>(*bmat % *grid_->gradz()->cut(offset, offset+b->nbasis()));

    double sum[3] = {0.0};
    for (size_t i = 0; i != grid_->size(); ++i) {
      for (int x = 0; x != 3; ++x)
        sum[x] += ddot_(mat->mdim(), d1mat[x]->element_ptr(0,i), 1, orb[0]->element_ptr(0,i), 1) * grid_->weight(i) * vxc[i];
    }

    if (!func->lda()) {
      array<shared_ptr<const Matrix>,6> d2mat;
      for (int i = 0; i != 6; ++i)
        d2mat[i] = make_shared<const Matrix>(*bmat % *grad2[i]->cut(offset, offset+b->nbasis()));

      unique_ptr<double[]> tmp2(new double[mat->mdim()]);
      for (size_t i = 0; i != grid_->size(); ++i) {
        // first term
        fill_n(tmp2.get(), mat->mdim(), 0.0);
        daxpy_(mat->mdim(), rhox[i], d2mat[0]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoy[i], d2mat[1]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoz[i], d2mat[3]->element_ptr(0,i), 1, tmp2.get(), 1);
        sum[0] += ddot_(mat->mdim(), tmp2.get(), 1, orb[0]->element_ptr(0,i), 1) * grid_->weight(i) * (2*vxc[i+grid_->size()]);
        fill_n(tmp2.get(), mat->mdim(), 0.0);
        daxpy_(mat->mdim(), rhox[i], d2mat[1]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoy[i], d2mat[2]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoz[i], d2mat[4]->element_ptr(0,i), 1, tmp2.get(), 1);
        sum[1] += ddot_(mat->mdim(), tmp2.get(), 1, orb[0]->element_ptr(0,i), 1) * grid_->weight(i) * (2*vxc[i+grid_->size()]);
        fill_n(tmp2.get(), mat->mdim(), 0.0);
        daxpy_(mat->mdim(), rhox[i], d2mat[3]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoy[i], d2mat[4]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoz[i], d2mat[5]->element_ptr(0,i), 1, tmp2.get(), 1);
        sum[2] += ddot_(mat->mdim(), tmp2.get(), 1, orb[0]->element_ptr(0,i), 1) * grid_->weight(i) * (2*vxc[i+grid_->size()]);
        // second term
        fill_n(tmp2.get(), mat->mdim(), 0.0);
        daxpy_(mat->mdim(), rhox[i], orb[1]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoy[i], orb[2]->element_ptr(0,i), 1, tmp2.get(), 1);
        daxpy_(mat->mdim(), rhoz[i], orb[3]->element_ptr(0,i), 1, tmp2.get(), 1);
        for (int x = 0; x != 3; ++x)
          sum[x] += ddot_(mat->mdim(), tmp2.get(), 1, d1mat[x]->element_ptr(0,i), 1) * grid_->weight(i) * (2*vxc[i+grid_->size()]);
      }
    }

    out->element(0, n) += -4.0*sum[0];
    out->element(1, n) += -4.0*sum[1];
    out->element(2, n) += -4.0*sum[2];

    offset += b->nbasis();
    ++n;
  }

  return out;
}


constexpr double a_stratmann__ = 0.64;

double DFTGrid_base::fuzzy_cell(shared_ptr<const Atom> atom, array<double,3>&& xyz) const {
  int fuzzy = -1;
  shared_ptr<StackMem> stack = resources__->get();
  double* const total = stack->get(mol_->natom());
  fill_n(total, mol_->natom(), 1.0);

  int i = 0;
  for (auto b = mol_->atoms().begin(); b != mol_->atoms().end(); ++b, ++i) {
    const double distbg = (*b)->distance(xyz);
    int j = i+1;
    for (auto c = b+1; c != mol_->atoms().end(); ++c, ++j) {
      const double distbc = (*b)->distance(*c);
      const double distcg = (*c)->distance(xyz);
      // Stratmann CPL 1996
      double nuij = (distbg - distcg) / distbc / a_stratmann__;
      const double nuij2 = pow(nuij,2);
      const double fac = (nuij >= -1 ? (nuij < 1 ?  0.5-0.5*(35.0/16.0*nuij*(1.0-nuij2*(1.0-21.0/35.0*nuij2*(1.0-5.0/21.0*nuij2)))) : 0.0) : 1.0);
      total[i] *= fac;
      total[j] *= 1.0-fac;
    }
    // TODO threshold is still hardwired (not a good practice)
    if ((*b)->distance(atom) < 1.0e-3) fuzzy = i;
  }

  if (fuzzy == -1)
    throw runtime_error("grid and atoms do not match with each other");

  const double out = total[fuzzy] / accumulate(total, total+mol_->natom(), 0.0); // Eq. 22
  stack->release(mol_->natom(), total);
  resources__->release(stack);
  return out;
}


namespace bagel {
class FuzzyTask {
  protected:
    shared_ptr<Matrix> data;
    shared_ptr<const Atom> atom;
    double xg;
    double yg;
    double zg;
    double coeff;
    DFTGrid_base* parent;
    const int n;
  public:
    FuzzyTask(shared_ptr<Matrix> d, shared_ptr<const Atom> a, double x, double y, double z, double c, DFTGrid_base* ptr, const int i)
     : data(d), atom(a), xg(x), yg(y), zg(z), coeff(c), parent(ptr), n(i) { }

    void compute() {
      const double weight = coeff * parent->fuzzy_cell(atom, array<double,3>{{xg, yg, zg}});
      data->element(0, n) = xg;
      data->element(1, n) = yg;
      data->element(2, n) = zg;
      data->element(3, n) = weight;
    }

};
}


void DFTGrid_base::add_grid(const int nrad, const int nang, const unique_ptr<double[]>& r_ch, const unique_ptr<double[]>& w_ch,
                            const unique_ptr<double[]>& x, const unique_ptr<double[]>& y, const unique_ptr<double[]>& z, const unique_ptr<double[]>& w) {

  const int ngrid = nrad*nang;
  const int nprev = grid_ ? grid_->size() : 0;

  auto combined = make_shared<Matrix>(4, nprev+ngrid*mol_->natom());
  if (nprev)
    copy_n(grid_->data()->data(), 4*nprev, combined->data());

  TaskQueue<FuzzyTask> tasks(mol_->natom()*nrad*nang);

  int cnt = nprev;
  for (auto& a : mol_->atoms()) {
    const double rbs = a->radius();

    double rib = 1.0e+10;
    for (auto& b : mol_->atoms())
      if (a != b)
        rib = min(rib, a->distance(b));

    for (int i = 0; i != nrad; ++i) {
      const double rr = r_ch[i] * rbs;
      if (rr < (0.5-0.5*a_stratmann__)*rib) {
        for (int j = 0; j != nang; ++j) {
          const double xg = x[j] * rr + a->position(0);
          const double yg = y[j] * rr + a->position(1);
          const double zg = z[j] * rr + a->position(2);
          combined->element(0, cnt) = xg;
          combined->element(1, cnt) = yg;
          combined->element(2, cnt) = zg;
          combined->element(3, cnt) = w[j]*w_ch[i]*pow(rbs,3)*4.0*pi__;
          ++cnt;
        }
      } else {
        for (int j = 0; j != nang; ++j) {
          const double xg = x[j] * rr + a->position(0);
          const double yg = y[j] * rr + a->position(1);
          const double zg = z[j] * rr + a->position(2);
          tasks.emplace_back(combined, a, xg, yg, zg, w[j]*w_ch[i]*pow(rbs,3)*4.0*pi__, this, cnt++);
        }
      }
    }
  }
  tasks.compute();

  shared_ptr<const Matrix> o = combined;
  grid_ = make_shared<Grid>(mol_, o);

}


void DFTGrid_base::remove_redgrid() {
  int size = 0;
  for (int i = 0; i != grid_->size(); ++i)
    if (grid_->data()->element(3, i) > grid_thresh_)
      ++size;
  auto out = make_shared<Matrix>(4, size);

  if (size < grid_->size())
    cout << "    * Removing " << grid_->size()-size << " points whose weight is below " << scientific << setprecision(2) << grid_thresh_ << endl << fixed;

  size = 0;
  for (int i = 0; i != grid_->size(); ++i)
    if (grid_->data()->element(3, i) > grid_thresh_)
      copy_n(grid_->data()->element_ptr(0, i), 4, out->element_ptr(0,size++));

  shared_ptr<const Matrix> o = out;
  grid_ = make_shared<Grid>(mol_, o);

  cout <<  "    * Grid points: " << size << endl << endl;
}


// grid without 'pruning'. Becke's original mapping
BLGrid::BLGrid(const size_t nrad, const size_t nang, shared_ptr<const Molecule> mol) : DFTGrid_base(mol) {
  // construct Lebedev grid
  unique_ptr<double[]> x(new double[nang]);
  unique_ptr<double[]> y(new double[nang]);
  unique_ptr<double[]> z(new double[nang]);
  unique_ptr<double[]> w(new double[nang]);
  lebedev.root(nang, x.get(), y.get(), z.get(), w.get());

  // construct Chebyshev grid
  unique_ptr<double[]> r_ch(new double[nrad]);
  unique_ptr<double[]> w_ch(new double[nrad]);
  for (int i = 0; i != nrad; ++i) {
    const double t = cos((i+1)*pi__/(nrad+1));
    r_ch[i] = (1.0+t)/(1.0-t);
    w_ch[i] = 2.0 / pow(1.0-t, 2.0)                  // due to mapping from [0,infty) to [-1, 1]
            * pi__/(nrad+1)*sin((i+1)*pi__/(nrad+1)) // Gauss-Chebyshev weight
            * r_ch[i]*r_ch[i];                       // due to r^2 in the spherical coordinate integration
  }

  add_grid(nrad, nang, r_ch, w_ch, x, y, z, w);
  remove_redgrid();
  grid_->init();
}


TALGrid::TALGrid(const size_t nrad, const size_t nang, shared_ptr<const Molecule> mol) : DFTGrid_base(mol) {
  // construct Lebedev grid
  unique_ptr<double[]> x(new double[nang]);
  unique_ptr<double[]> y(new double[nang]);
  unique_ptr<double[]> z(new double[nang]);
  unique_ptr<double[]> w(new double[nang]);
  lebedev.root(nang, x.get(), y.get(), z.get(), w.get());

  // construct Chebyshev grid
  unique_ptr<double[]> r_ch(new double[nrad]);
  unique_ptr<double[]> w_ch(new double[nrad]);
  for (int i = 0; i != nrad; ++i) {
    const double t = cos((i+1)*pi__/(nrad+1));
    r_ch[i] = 1.0/log(2.0)*pow(1.0+t, 0.6)*log(2.0/(1.0-t));
    w_ch[i] = (0.6*r_ch[i]/(1.0+t) + 1.0/log(2.0)*pow(1.0+t,0.6)/(1-t)) // due to mapping from [0,infty) to [-1, 1]
            * pi__/(nrad+1)*sin((i+1)*pi__/(nrad+1)) // Gauss-Chebyshev weight
            * r_ch[i]*r_ch[i];                       // due to r^2 in the spherical coordinate integration
  }

  add_grid(nrad, nang, r_ch, w_ch, x, y, z, w);
  remove_redgrid();
  grid_->init();
}


DefaultGrid::DefaultGrid(shared_ptr<const Molecule> mol) : DFTGrid_base(mol) {
  // the default radial grid has 75 points
  const int nrad = 75;
  // construct Chebyshev grid
  unique_ptr<double[]> r_ch(new double[nrad]);
  unique_ptr<double[]> w_ch(new double[nrad]);
  for (int i = 0; i != nrad; ++i) {
    const double t = cos((i+1)*pi__/(nrad+1));
    r_ch[i] = 1.0/log(2.0)*pow(1.0+t, 0.6)*log(2.0/(1.0-t));
    w_ch[i] = (0.6*r_ch[i]/(1.0+t) + 1.0/log(2.0)*pow(1.0+t,0.6)/(1-t)) // due to mapping from [0,infty) to [-1, 1]
            * pi__/(nrad+1)*sin((i+1)*pi__/(nrad+1)) // Gauss-Chebyshev weight
            * r_ch[i]*r_ch[i];                       // due to r^2 in the spherical coordinate integration
  }

  // start, fence, nang
  vector<tuple<int, int, int>> map;
  // TODO Decide how to partition
  map.push_back(make_tuple(0, 8, 194));
  map.push_back(make_tuple(8, 45, 302));
  map.push_back(make_tuple(45, 50, 194));
  map.push_back(make_tuple(50, 55, 110));
  map.push_back(make_tuple(55, 60, 50));
  map.push_back(make_tuple(60, 70, 38));
  map.push_back(make_tuple(70, 75, 6));
  for (auto& i : map) {
    const int nang = get<2>(i);
    unique_ptr<double[]> rr(new double[get<1>(i)-get<0>(i)]);
    unique_ptr<double[]> ww(new double[get<1>(i)-get<0>(i)]);
    copy(r_ch.get()+get<0>(i), r_ch.get()+get<1>(i), rr.get());
    copy(w_ch.get()+get<0>(i), w_ch.get()+get<1>(i), ww.get());
    // construct Lebedev grid
    unique_ptr<double[]> x(new double[nang]);
    unique_ptr<double[]> y(new double[nang]);
    unique_ptr<double[]> z(new double[nang]);
    unique_ptr<double[]> w(new double[nang]);
    lebedev.root(nang, x.get(), y.get(), z.get(), w.get());
    add_grid(get<1>(i)-get<0>(i), nang, rr, ww, x, y, z, w);
  }

  remove_redgrid();
  grid_->init();
}
