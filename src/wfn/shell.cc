//
// BAGEL - Parallel electron correlation program.
// Filename: shell.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#include <src/osint/overlapbatch.h>
#include <src/osint/momentbatch.h>
#include <src/rysint/carsphlist.h>
#include <src/wfn/shell.h>
#include <src/rysint/intparam.h> // for Cartesian to spherical transformation
#include <src/rysint/carsphlist.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace bagel;

static const CarSphList carsphlist;

Shell::Shell(const bool sph, const array<double,3>& _position, int _ang, const vector<double>& _expo,
                       const vector<vector<double>>& _contr,  const vector<pair<int, int>>& _range)
 : spherical_(sph), position_(_position), angular_number_(_ang),
   exponents_(_expo), contractions_(_contr), contraction_ranges_(_range), dummy_(false), relativistic_(false) {

  contraction_lower_.reserve(_range.size());
  contraction_upper_.reserve(_range.size());
  for (auto piter = _range.begin(); piter != _range.end(); ++piter) {
    contraction_lower_.push_back(piter->first);
    contraction_upper_.push_back(piter->second);
  }

  if (spherical_)
    nbasis_ = (angular_number_*2+1) * num_contracted();
  else
    nbasis_ = (angular_number_+1) * (angular_number_+2) / 2 * num_contracted();

}


Shell::Shell(const bool sph) : spherical_(sph), position_{{0.0,0.0,0.0}}, angular_number_(0), exponents_(1,0.0), contraction_ranges_(1,make_pair(0,1)),
        dummy_(true) {
  contractions_.push_back(vector<double>(1,1.0));
  contraction_lower_.push_back(0);
  contraction_upper_.push_back(1);
}


Shell::~Shell() {

}


std::shared_ptr<const Shell> Shell::move_atom(const array<double,3>& displacement) const {
  std::shared_ptr<Shell> out(new Shell(*this));
  out->position_[0] += displacement[0];
  out->position_[1] += displacement[1];
  out->position_[2] += displacement[2];
  return out;
}


std::shared_ptr<const Shell> Shell::move_atom(const double* displacement) const {
  std::shared_ptr<Shell> out(new Shell(*this));
  out->position_[0] += displacement[0];
  out->position_[1] += displacement[1];
  out->position_[2] += displacement[2];
  return out;
}


const string Shell::show() const {
  stringstream ss;
  ss << "position: ";
  ss << position_[0] << " " << position_[1] << " "  << position_[2] << endl;
  ss << "angular: "  << angular_number_ << endl;
  ss << "exponents: ";
  for (int i = 0; i != exponents_.size(); ++i) {
    ss << " " << exponents_[i];
  }
  ss << endl;
  for (int i = 0; i != contractions_.size(); ++i) {
    ss << " (" << contraction_ranges_[i].first << "," << contraction_ranges_[i].second << ") ";
    for (int j = contraction_ranges_[i].first; j != contraction_ranges_[i].second; ++j)
      ss << contractions_[i][j] << " ";
  }

  return ss.str();
}


bool Shell::operator==(const Shell& o) const {
  bool out = true;
  out &= spherical_ == o.spherical_;
  out &= position_ == o.position_;
  out &= angular_number_ == o.angular_number_;
  out &= exponents_ == o.exponents_;
  out &= contractions_ == o.contractions_;
  out &= contraction_ranges_ == o.contraction_ranges_;
  out &= dummy_ == o.dummy_;
  out &= contraction_upper_ == o.contraction_upper_;
  out &= contraction_lower_ == o.contraction_lower_;
  out &= nbasis_ == o.nbasis_;
  return out;
}


// returns uncontracted cartesian shell with one higher or lower angular number if inc is + or - 1 respectively
shared_ptr<const Shell> Shell::kinetic_balance_uncont(int inc) const {
  assert(abs(inc)==1);
  int i = 0;
  vector<vector<double>> conts;
  vector<pair<int, int>> ranges;
  for (auto e = exponents_.begin(); e != exponents_.end(); ++e, ++i) {
    vector<double> cont(exponents_.size(), 0);
    cont[i] = 1.0;
    conts.push_back(cont);
    ranges.push_back(make_pair(i,i+1));
  }
  return angular_number_+inc < 0 ? shared_ptr<const Shell>() : shared_ptr<const Shell>(new Shell(false, position_, angular_number_+inc, exponents_, conts, ranges));
}


shared_ptr<const Shell> Shell::cartesian_shell() const {
  shared_ptr<Shell> out(new Shell(false, position_, angular_number_, exponents_, contractions_, contraction_ranges_));
  return out;
}


void Shell::init_relativistic() {
  relativistic_ = true;
  aux_dec_ = kinetic_balance_uncont(-1);
  aux_inc_ = kinetic_balance_uncont(1);

  // overlap = S^-1 between auxiliary functions
  shared_ptr<const Matrix> overlap = overlap_compute_();

  // small is a transformation matrix (x,y,z components)
  small_ = moment_compute_(overlap);
}


shared_ptr<const Matrix> Shell::overlap_compute_() const {

  const int asize_inc = aux_inc_->nbasis();
  const int asize_dec = aux_dec_ ? aux_dec_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  shared_ptr<Matrix> overlap(new Matrix(a,a, true));

  {
    OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_, aux_inc_}});
    ovl.compute();
    for (int i = 0; i != aux_inc_->nbasis(); ++i)
      copy_n(ovl.data() + i*asize_inc, asize_inc, overlap->element_ptr(0,i));
  }
  if (aux_dec_) {
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_dec_, aux_dec_}});
      ovl.compute();
      for (int i = 0; i != asize_dec; ++i) 
        copy_n(ovl.data() + i*asize_dec, asize_dec, overlap->element_ptr(asize_inc, i+asize_inc));
    }
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_, aux_dec_}});
      ovl.compute();
      for (int i = 0; i != asize_dec; ++i)
        for (int j = 0; j != asize_inc; ++j)
          overlap->element(j,i+asize_inc) = overlap->element(i+asize_inc,j) = *(ovl.data()+j+asize_inc*i);
    }
  }

  return overlap;
}


array<shared_ptr<const Matrix>,3> Shell::moment_compute_(const shared_ptr<const Matrix> overlap) const {
  const int ssize = nbasis();
  const int asize_inc = aux_inc_->nbasis();
  const int asize_dec = aux_dec_ ? aux_dec_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  shared_ptr<MomentBatch> coeff0(new MomentBatch(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_inc_}}));
  coeff0->compute();

  shared_ptr<MomentBatch> coeff1;
  if (aux_dec_) {
    coeff1 = shared_ptr<MomentBatch>(new MomentBatch(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_dec_}}));
    coeff1->compute();
  } else {
    // just to run. coeff1 is not referenced in the code
    coeff1 = coeff0;
  }

  const double* carea0 = coeff0->data();
  const double* carea1 = coeff1->data();

  shared_ptr<Matrix> tmparea(new Matrix(ssize,a, true));
  array<shared_ptr<const Matrix>,3> out;

  const static CarSphList carsphlist;
  for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
    if (spherical_) {
      const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
      const int nloop = num_contracted() * asize_inc; 
      carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea->data());
    } else {
      assert(coeff0->size_block() == asize_inc*ssize);
      copy(carea0, carea0+coeff0->size_block(), tmparea->data());
    }
    if (aux_dec_) {
      if (spherical_) {
        const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
        const int nloop = num_contracted() * asize_dec; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea->data()+asize_inc*ssize);
      } else {
        assert(coeff1->size_block() == asize_dec*ssize);
        copy(carea1, carea1+coeff1->size_block(), tmparea->data()+asize_inc*ssize);
        }
    }

    out[i] = tmparea->transpose()->solve(overlap, overlap->ndim());
  }
  return out;
}


// In DFT we want to compute values of basis functions on grid
void Shell::compute_grid_value(double* b, double* dx, double* dy, double* dz, const double& x, const double& y, const double& z) const {
  const double rr = x*x+y*y+z*z;
  auto range = contraction_ranges_.begin();
  double tmp0[50];
  double tmpx[50];
  double tmpy[50];
  double tmpz[50];
  assert(50 > ANG_HRR_END*ANG_HRR_END);

  const int nxyz = nbasis() / num_contracted();
  const int nang = angular_number();
  const int index = nang * ANG_HRR_END;

  for (auto& i : contractions_) {
    double exp0 = 0.0;
    double expx = 0.0;
    double expy = 0.0;
    double expz = 0.0;
    for (int j = range->first; j != range->second; ++j) { 
      const double tmp = i[j]*exp(-exponents_[j]*rr);
      exp0 += tmp; 
      expx += -2.0*exponents_[j]*x*tmp;
      expy += -2.0*exponents_[j]*y*tmp;
      expz += -2.0*exponents_[j]*z*tmp;
    }
    for (int iz = 0, ixyz = 0; iz <= nang; ++iz) {
      for (int iy = 0; iy <= nang - iz; ++iy, ++ixyz) {
        const int ix = nang - iy - iz;
        const double cart = pow(x,ix)*pow(y,iy)*pow(z,iz);
        tmp0[ixyz] = cart*exp0;
        tmpx[ixyz] = (ix != 0 ? ix*pow(x,ix-1)*pow(y,iy)*pow(z,iz)*exp0 : 0.0) + cart*expx;
        tmpy[ixyz] = (iy != 0 ? iy*pow(x,ix)*pow(y,iy-1)*pow(z,iz)*exp0 : 0.0) + cart*expy;
        tmpz[ixyz] = (iz != 0 ? iz*pow(x,ix)*pow(y,iy)*pow(z,iz-1)*exp0 : 0.0) + cart*expz;
      }
    }
    if (spherical_) {
      carsphlist.carsphfunc_call(index, 1, tmp0, b);
      carsphlist.carsphfunc_call(index, 1, tmpx, dx);
      carsphlist.carsphfunc_call(index, 1, tmpy, dy);
      carsphlist.carsphfunc_call(index, 1, tmpz, dz);
    } else {
      copy_n(tmp0, nxyz, b);
      copy_n(tmpx, nxyz, dx);
      copy_n(tmpy, nxyz, dy);
      copy_n(tmpz, nxyz, dz);
    }
    b  += nxyz;
    dx += nxyz;
    dy += nxyz;
    dz += nxyz;
    ++range;
  }
}
