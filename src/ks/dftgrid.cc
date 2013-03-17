//
// BAGEL - Parallel electron correlation program.
// Filename: dftgrid.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <numeric>
#include <src/util/f77.h>
#include <src/ks/dftgrid.h>
#include <src/ks/lebedevlist.h>
#include <src/util/constants.h>
#include <src/parallel/mpi_interface.h>
#include <src/ks/xcfunc.h>

using namespace std;
using namespace bagel;

const static LebedevList lebedev;


void DFTGridPoint::init() {
  const int ngrid = size();
  basis_ = shared_ptr<Matrix>(new Matrix(geom_->nbasis(), ngrid));
  gradx_ = shared_ptr<Matrix>(new Matrix(geom_->nbasis(), ngrid));
  grady_ = shared_ptr<Matrix>(new Matrix(geom_->nbasis(), ngrid));
  gradz_ = shared_ptr<Matrix>(new Matrix(geom_->nbasis(), ngrid));

  for (size_t g = 0; g != ngrid; ++g) {
    int pos = 0;
    for (auto& i : geom_->atoms()) {
      // xyz coordinate relative to the atom i
      const double x = data_->element(0,g) - i->position(0);
      const double y = data_->element(1,g) - i->position(1);
      const double z = data_->element(2,g) - i->position(2);
      for (auto& j : i->shells()) {
        // angular number
        j->compute_grid_value(basis_->element_ptr(pos, g), gradx_->element_ptr(pos, g), grady_->element_ptr(pos, g), gradz_->element_ptr(pos, g), x, y, z);
        pos += j->nbasis();
      }
    }
  }
}


shared_ptr<const Matrix> DFTGrid_base::compute_rho_sigma(shared_ptr<const XCFunc> func, shared_ptr<const Matrix> mat,
                                                         unique_ptr<double[]>& rho, unique_ptr<double[]>& sigma,
                                                         unique_ptr<double[]>& rhox, unique_ptr<double[]>& rhoy, unique_ptr<double[]>& rhoz) const {
  shared_ptr<Matrix> orb(new Matrix(*mat % *grid_->basis()));
  if (func->lda()) { 
    assert(orb->mdim() == grid_->size());
    for (size_t i = 0; i != orb->mdim(); ++i) {
      rho[i] = 2*ddot_(orb->ndim(), orb->element_ptr(0, i), 1, orb->element_ptr(0, i), 1);
    }
  } else {
    shared_ptr<Matrix> orb(new Matrix(*mat % *grid_->basis()));
    shared_ptr<Matrix> orbx(new Matrix(*mat % *grid_->gradx()));
    shared_ptr<Matrix> orby(new Matrix(*mat % *grid_->grady()));
    shared_ptr<Matrix> orbz(new Matrix(*mat % *grid_->gradz()));
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
  }
  return orb;
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

  unique_ptr<double[]> exc;
  unique_ptr<double[]> vxc;
  tie(exc, vxc) = func->compute_exc_vxc(grid_->size(), rho, sigma); 
  time.tick_print("exc");

  shared_ptr<Matrix> out(new Matrix(geom_->nbasis(), geom_->nbasis()));
  double en = 0.0;
  {
    shared_ptr<Matrix> scal(new Matrix(*grid_->basis()));
    for (size_t i = 0; i != scal->mdim(); ++i) {
      dscal_(scal->ndim(), vxc[i]*grid_->weight(i), scal->element_ptr(0, i), 1);
      en += exc[i] * rho[i] * grid_->weight(i);
    }
    *out += *grid_->basis() ^ *scal;
  }
  if (!func->lda()) {
    {
      shared_ptr<Matrix> scal(new Matrix(*grid_->basis()));
      for (size_t i = 0; i != scal->mdim(); ++i)
        dscal_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhox[i], scal->element_ptr(0, i), 1);
      *out += *grid_->gradx() ^ *scal;
    }
    {
      shared_ptr<Matrix> scal(new Matrix(*grid_->basis()));
      for (size_t i = 0; i != scal->mdim(); ++i)
        dscal_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhoy[i], scal->element_ptr(0, i), 1);
      *out += *grid_->grady() ^ *scal;
    }
    {
      shared_ptr<Matrix> scal(new Matrix(*grid_->basis()));
      for (size_t i = 0; i != scal->mdim(); ++i)
        dscal_(scal->ndim(), 4*vxc[i+grid_->size()]*grid_->weight(i)*rhoz[i], scal->element_ptr(0, i), 1);
      *out += *grid_->gradz() ^ *scal;
    }
  }
  out->symmetrize();

  time.tick_print("contraction");
  return make_tuple(out, en);
}


shared_ptr<const GradFile> DFTGrid_base::compute_xcgrad(shared_ptr<const XCFunc> func, shared_ptr<const Matrix> mat) const {
  shared_ptr<GradFile> out(new GradFile(geom_->natom()));

  unique_ptr<double[]> rho(new double[grid_->size()]);
  unique_ptr<double[]> sigma, rhox, rhoy, rhoz;
  if (!func->lda()) {
    sigma = unique_ptr<double[]>(new double[grid_->size()]);
    rhox  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoy  = unique_ptr<double[]>(new double[grid_->size()]);
    rhoz  = unique_ptr<double[]>(new double[grid_->size()]);
  }

  shared_ptr<const Matrix> orb = compute_rho_sigma(func, mat, rho, sigma, rhox, rhoy, rhoz);

  unique_ptr<double[]> vxc = func->compute_vxc(grid_->size(), rho, sigma);

  // loop over target atom
  size_t offset = 0;
  int n = 0;
  for (auto& b : geom_->atoms()) {
    shared_ptr<const Matrix> bmat = mat->cut(offset, offset+b->nbasis());
    shared_ptr<const Matrix> xmat(new Matrix(*bmat % *grid_->gradx()->cut(offset, offset+b->nbasis())));
    shared_ptr<const Matrix> ymat(new Matrix(*bmat % *grid_->grady()->cut(offset, offset+b->nbasis())));
    shared_ptr<const Matrix> zmat(new Matrix(*bmat % *grid_->gradz()->cut(offset, offset+b->nbasis())));

    double sum[3] = {0.0};
    for (size_t i = 0; i != grid_->size(); ++i) {
      sum[0] += ddot_(mat->mdim(), xmat->element_ptr(0,i), 1, orb->element_ptr(0,i), 1) * grid_->weight(i) * vxc[i];
      sum[1] += ddot_(mat->mdim(), ymat->element_ptr(0,i), 1, orb->element_ptr(0,i), 1) * grid_->weight(i) * vxc[i];
      sum[2] += ddot_(mat->mdim(), zmat->element_ptr(0,i), 1, orb->element_ptr(0,i), 1) * grid_->weight(i) * vxc[i];
    }
    out->data(0, n) += -4.0*sum[0];
    out->data(1, n) += -4.0*sum[1];
    out->data(2, n) += -4.0*sum[2];

    offset += b->nbasis();
    ++n;
  }

  return out; 
}


double DFTGrid_base::fuzzy_cell(shared_ptr<const Atom> atom, array<double,3>&& xyz) const {
  double fuzzy = -1.0;
  double total = 0.0;
  for (auto& b : geom_->atoms()) {
    const double rbs1 = b->radius();
    double tmp = 1.0;
    for (auto& c : geom_->atoms()) {
      if (b != c) {
        const double rbs2 = c->radius();
        const double xi = sqrt(rbs1/rbs2); // sqrt. see JCP 102, 346 (1995)
        const double uij = (xi-1.0)/(xi+1.0);
        const double aij = uij / (uij*uij-1.0);

        const double distbc = b->distance(c); 
        const double distbg = b->distance(xyz);
        const double distcg = c->distance(xyz);
        const double muij = (distbg - distcg) / distbc;

        // see Becke's appendix
        double nuij = muij + aij*(1.0-muij*muij); // eq. a2
        for (int i = 0; i != 3; ++i)
          nuij = (1.5-0.5*nuij*nuij)*nuij; // eq. 19

        tmp *= 0.5*(1.0-nuij); // eq. 21
      }
    }
    // TODO threshold is still hardwired (not a good practice)
    if (b->distance(atom) < 1.0e-3) fuzzy = tmp; 
    total += tmp;
  }

  if (fuzzy == -1.0)
    throw runtime_error("grid and atoms do not match with each other");

  return fuzzy / total; // Eq. 22
}


void DFTGrid_base::add_grid(const int nrad, const int nang, const unique_ptr<double[]>& r_ch, const unique_ptr<double[]>& w_ch,
                            const unique_ptr<double[]>& x, const unique_ptr<double[]>& y, const unique_ptr<double[]>& z, const unique_ptr<double[]>& w) {

  const int ngrid = nrad*nang;
  const int nprev = grid_ ? grid_->size() : 0;
  shared_ptr<Matrix> data(new Matrix(4, nprev+ngrid*geom_->natom()));
  if (nprev)
    copy_n(grid_->data()->data(), 4*nprev, data->data());

  int cnt = nprev;
  for (auto& a : geom_->atoms()) {
    const double rbs = a->radius();
    for (int i = 0; i != nrad; ++i) {
      for (int j = 0; j != nang; ++j) {
        const double xg = x[j] * r_ch[i] * rbs + a->position(0);
        const double yg = y[j] * r_ch[i] * rbs + a->position(1);
        const double zg = z[j] * r_ch[i] * rbs + a->position(2);
        double weight = w[j] * w_ch[i] * pow(rbs,3) * 4.0*pi__ * fuzzy_cell(a, array<double,3>{{xg, yg, zg}});

        // set to data 
        if (weight > grid_thresh_/(nang*nrad)) {
          data->element(0, cnt) = xg;
          data->element(1, cnt) = yg;
          data->element(2, cnt) = zg;
          data->element(3, cnt) = weight;
          ++cnt;
        }
      }
    }
  }
  shared_ptr<const Matrix> o = data->slice(0,cnt);
  grid_ = shared_ptr<const DFTGridPoint>(new DFTGridPoint(geom_, o));

}


// grid without 'pruning'. Becke's original mapping
BLGrid::BLGrid(const size_t nrad, const size_t nang, shared_ptr<const Geometry> geom) : DFTGrid_base(geom) {
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
}


TALGrid::TALGrid(const size_t nrad, const size_t nang, shared_ptr<const Geometry> geom) : DFTGrid_base(geom) {
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
}


DefaultGrid::DefaultGrid(shared_ptr<const Geometry> geom) : DFTGrid_base(geom) {
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
  map.push_back(make_tuple(0,  10, 302));
  map.push_back(make_tuple(10, nrad, 302));
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
}
