//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: shell.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/molecule/shell.h>
#include <src/molecule/moment_compute.h>
#include <src/util/math/matop.h>
#include <src/integral/carsphlist.h>

using namespace std;
using namespace bagel;

static const CarSphList carsphlist;

Shell::Shell(const bool sph, const array<double,3>& _position, int _ang, const vector<double>& _expo,
                       const vector<vector<double>>& _contr,  const vector<pair<int, int>>& _range)
 : Shell_base(sph, _position, _ang), exponents_(_expo), contractions_(_contr), contraction_ranges_(_range),
   dummy_(false), relativistic_(false), magnetism_(false), london_(false), vector_potential_{{0.0, 0.0, 0.0}}, magnetic_field_{{0.0, 0.0, 0.0}} {

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


Shell::Shell(const bool sph) : Shell_base(sph), exponents_{0.0}, contractions_{{1.0}}, contraction_ranges_{{0,1}},
                               dummy_(true), magnetism_(false), london_(false), vector_potential_{{0.0, 0.0, 0.0}}, magnetic_field_{{0.0, 0.0, 0.0}} {
  contraction_lower_.push_back(0);
  contraction_upper_.push_back(1);
}


shared_ptr<const Shell> Shell::move_atom(const array<double,3>& displacement) const {
  auto out = make_shared<Shell>(*this);
  out->position_[0] += displacement[0];
  out->position_[1] += displacement[1];
  out->position_[2] += displacement[2];
  return out;
}


shared_ptr<const Shell> Shell::move_atom(const double* displacement) const {
  auto out = make_shared<Shell>(*this);
  out->position_[0] += displacement[0];
  out->position_[1] += displacement[1];
  out->position_[2] += displacement[2];
  return out;
}

string Shell::show() const {
  stringstream ss;
  ss << "position: ";
  ss << position_[0] << " " << position_[1] << " "  << position_[2] << endl;
  ss << "angular: "  << angular_number_ << endl;
  if (magnetism_)
    ss << "vector potential: " << vector_potential_[0] << " " << vector_potential_[1] << " " << vector_potential_[2] << endl;
  ss << "exponents: ";
  for (int i = 0; i != exponents_.size(); ++i) {
    ss << " " << exponents_[i];
  }
  ss << endl;
  ss << "contraction coefficients: ";
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

vector<shared_ptr<const Shell>> Shell::split_if_possible(const size_t batchsize) const {
  vector<shared_ptr<const Shell>> out;
  // first see if there are disconnected shells
  const int nb = nbasis_ / contraction_upper_.size();
  assert(nbasis_%contraction_upper_.size() == 0);

  int smallest = 0;
  int largest = contraction_upper_.front();
  auto upper = contraction_upper_.begin();
  auto lower = contraction_lower_.begin();
  int nstart = 0;
  int nend = 0;
  while (1) {
    ++upper;
    ++lower;
    ++nend;
    // if this condition is met, we make a shell object
    if (upper == contraction_upper_.end() || (*lower >= largest && (nend-nstart)*nb >= batchsize)) {
      vector<double> expo(exponents_.begin()+smallest, exponents_.begin()+largest);
      vector<vector<double>> contr;
      vector<pair<int,int>>  range;
      for (int i = nstart; i != nend; ++i) {
        contr.push_back(vector<double>(contractions_[i].begin()+smallest, contractions_[i].end()));
        range.push_back({contraction_ranges_[i].first-smallest, contraction_ranges_[i].second-smallest});
      }
      out.push_back(make_shared<Shell>(spherical_, position_, angular_number_, expo, contr, range));
      smallest = *lower;
      nstart = nend;
      if (upper == contraction_upper_.end()) break;
    } else {
      if (*lower < smallest) throw runtime_error("This type of basis set is not supported yet.");
    }
    largest = max(*upper, largest);
  }
  assert(nend == contraction_upper_.size());
  return out;
}


// returns uncontracted cartesian shell with one higher or lower angular number if increment is + or - 1 respectively
template<int increment>
shared_ptr<const Shell> Shell::kinetic_balance_uncont() const {
  static_assert(increment==1||increment==0||increment==-1, "illegal call of Shell::kinetic_balance_uncont");
  int i = 0;
  vector<vector<double>> conts;
  vector<pair<int, int>> ranges;
  for (auto e = exponents_.begin(); e != exponents_.end(); ++e, ++i) {
    vector<double> cont(exponents_.size(), 0);
    cont[i] = 1.0;
    conts.push_back(cont);
    ranges.push_back({i,i+1});
  }
  auto out = angular_number_+increment < 0 ? nullptr : make_shared<Shell>(false, position_, angular_number_+increment, exponents_, conts, ranges);
  if (magnetism_ && angular_number_+increment >= 0) out->add_phase(vector_potential_, magnetic_field_, london_);
  return out;
}


void Shell::add_phase(const array<double,3>& phase_input, const array<double,3>& magnetic_field, const bool london) {
  assert(london || (phase_input[0] == 0.0 && phase_input[1] == 0.0 && phase_input[2] == 0.0));

  magnetism_ = true;
  london_ = london;
  vector_potential_ = phase_input;
  magnetic_field_ = magnetic_field;
}


shared_ptr<const Shell> Shell::cartesian_shell() const {
  auto out = make_shared<Shell>(false, position_, angular_number_, exponents_, contractions_, contraction_ranges_);
  if (magnetism_) out->add_phase(vector_potential_, magnetic_field_, london_);
  return out;
}


void Shell::init_relativistic() {
#ifndef COMPILE_J_ORB
  if (angular_number_ == 6) throw runtime_error("Relativistic calculations with i-type orbital basis functions require j-type integrals for the small component.  Recompile with -DCOMPILE_J_ORB to use this feature.");
#endif
  if (angular_number_ == 7) throw runtime_error("Relativistic codes cannot use j-type main basis functions, since k-type would be needed for the small component.");
  relativistic_ = true;
  aux_decrement_ = kinetic_balance_uncont<-1>();
  aux_increment_ = kinetic_balance_uncont<1>();

  // small is a transformation matrix (x,y,z components)
  small_ = MomentCompute::call(*this);
}


void Shell::init_relativistic(const array<double,3> magnetic_field, bool london) {
  assert(magnetism_);
#ifndef COMPILE_J_ORB
  if (angular_number_ == 6) throw runtime_error("Relativistic calculations with i-type orbital basis functions require j-type integrals for the small component.  Recompile with -DCOMPILE_J_ORB to use this feature.");
#endif
  if (angular_number_ == 7) throw runtime_error("Relativistic codes cannot use j-type main basis functions, since k-type would be needed for the small component.");
  relativistic_ = true;
  aux_decrement_ = kinetic_balance_uncont<-1>();
  aux_increment_ = kinetic_balance_uncont<1>();
  aux_same_ = london ? nullptr : kinetic_balance_uncont<0>();

  // zsmall is a transformation matrix (x,y,z components)
  zsmall_ = MomentCompute::call(*this, magnetic_field, london);
  for (int i = 0; i != 3; i++) zsmallc_[i] = zsmall_[i]->get_conjg();
}


// In DFT we want to compute values of basis functions on grid
void Shell::compute_grid_value(double* b, double* dx, double* dy, double* dz, const double& x, const double& y, const double& z) const {
  const bool dogradient = dx != nullptr && dy != nullptr && dz != nullptr;
  assert(dogradient || (dx == nullptr && dy == nullptr && dz == nullptr));

  const double rr = x*x+y*y+z*z;
  auto range = contraction_ranges_.begin();
  const int nang = angular_number();

  double tmp0[82];
  double tmpx[82];
  double tmpy[82];
  double tmpz[82];
  static_assert(82 > (ANG_HRR_END+1)*(ANG_HRR_END+1), "ANG_HRR_END is assumed to be 8");
  double powx[11], powy[11], powz[11];
  powx[0] =  powy[0] = powz[0] = 0.0;
  for (int i = 0; i != angular_number()+1; ++i) {
    powx[i+1] = pow(x, i);
    powy[i+1] = pow(y, i);
    powz[i+1] = pow(z, i);
  }

  const int nxyz = nbasis() / num_contracted();
  const int index = nang * ANG_HRR_END;

  for (auto& i : contractions_) {
    double exp0 = 0.0;
    double exp1 = 0.0;
    for (int j = range->first; j != range->second; ++j) {
      const double tmp = i[j]*exp(-exponents_[j]*rr);
      exp0 += tmp;
      exp1 -= 2.0*exponents_[j]*tmp;
    }
    const double expx = exp1*x;
    const double expy = exp1*y;
    const double expz = exp1*z;
    // TODO threshold hardwired
    const bool nonzero = fabs(exp0) >= 1.0e-14;
    if (nonzero && dogradient) {
      for (int iz = 0, ixyz = 0; iz <= nang; ++iz) {
        for (int iy = 0; iy <= nang - iz; ++iy, ++ixyz) {
          const int ix = nang - iy - iz;
          const double cart = powx[ix+1]*powy[iy+1]*powz[iz+1];
          tmp0[ixyz] = cart*exp0;
          tmpx[ixyz] = ix*powx[ix]*powy[iy+1]*powz[iz+1]*exp0 + cart*expx;
          tmpy[ixyz] = iy*powx[ix+1]*powy[iy]*powz[iz+1]*exp0 + cart*expy;
          tmpz[ixyz] = iz*powx[ix+1]*powy[iy+1]*powz[iz]*exp0 + cart*expz;
        }
      }
      if (spherical_ && index) {
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
    } else if (nonzero) {
      for (int iz = 0, ixyz = 0; iz <= nang; ++iz)
        for (int iy = 0; iy <= nang - iz; ++iy, ++ixyz) {
          const int ix = nang - iy - iz;
          const double cart = powx[ix+1]*powy[iy+1]*powz[iz+1];
          tmp0[ixyz] = cart*exp0;
        }
      if (spherical_ && index)
        carsphlist.carsphfunc_call(index, 1, tmp0, b);
      else
        copy_n(tmp0, nxyz, b);
    }
    b  += nxyz;
    if (dogradient) {
      dx += nxyz;
      dy += nxyz;
      dz += nxyz;
    }
    ++range;
  }
}


// In DFT we want to compute values of basis functions on grid
void Shell::compute_grid_value_deriv2(double* bxx, double* bxy, double* byy, double* bxz, double* byz, double* bzz,
                                      const double& x, const double& y, const double& z) const {
  const double rr = x*x+y*y+z*z;
  auto range = contraction_ranges_.begin();
  double tmp[6][65];
  assert(65 > ANG_HRR_END*ANG_HRR_END);

  const int nxyz = nbasis() / num_contracted();
  const int nang = angular_number();
  const int index = nang * ANG_HRR_END;

  for (auto& i : contractions_) {
    double exp0 = 0.0;
    double expx = 0.0;
    double expy = 0.0;
    double expz = 0.0;
    double expxx = 0.0;
    double expxy = 0.0;
    double expyy = 0.0;
    double expxz = 0.0;
    double expyz = 0.0;
    double expzz = 0.0;
    for (int j = range->first; j != range->second; ++j) {
      const double tmp = i[j]*exp(-exponents_[j]*rr);
      exp0 += tmp;
      expx += -2.0*exponents_[j]*x*tmp;
      expy += -2.0*exponents_[j]*y*tmp;
      expz += -2.0*exponents_[j]*z*tmp;
      expxx += -2.0*exponents_[j]*tmp + pow(2.0*exponents_[j]*x,2)*tmp;
      expxy += 4.0*pow(exponents_[j],2)*x*y*tmp;
      expyy += -2.0*exponents_[j]*tmp + pow(2.0*exponents_[j]*y,2)*tmp;
      expxz += 4.0*pow(exponents_[j],2)*x*z*tmp;
      expyz += 4.0*pow(exponents_[j],2)*y*z*tmp;
      expzz += -2.0*exponents_[j]*tmp + pow(2.0*exponents_[j]*z,2)*tmp;
    }
    for (int iz = 0, ixyz = 0; iz <= nang; ++iz) {
      for (int iy = 0; iy <= nang - iz; ++iy, ++ixyz) {
        const int ix = nang - iy - iz;
        const double cart = pow(x,ix)*pow(y,iy)*pow(z,iz);
        tmp[0][ixyz] = (ix > 1 ? ix*(ix-1)*pow(x,ix-2)*pow(y,iy)*pow(z,iz)*exp0 : 0)
                     + (ix > 0 ? 2.0*ix*pow(x,ix-1)*pow(y,iy)*pow(z,iz)*expx : 0)
                     + cart*expxx;
        tmp[1][ixyz] = (ix > 0 && iy > 0 ? ix*iy*pow(x,ix-1)*pow(y,iy-1)*pow(z,iz)*exp0 : 0)
                     + (ix > 0 ? ix*pow(x,ix-1)*pow(y,iy)*pow(z,iz)*expy : 0)
                     + (iy > 0 ? iy*pow(x,ix)*pow(y,iy-1)*pow(z,iz)*expx : 0)
                     + cart*expxy;
        tmp[2][ixyz] = (iy > 1 ? iy*(iy-1)*pow(x,ix)*pow(y,iy-2)*pow(z,iz)*exp0 : 0)
                     + (iy > 0 ? 2.0*iy*pow(x,ix)*pow(y,iy-1)*pow(z,iz)*expy : 0)
                     + cart*expyy;
        tmp[3][ixyz] = (ix > 0 && iz > 0 ? ix*iz*pow(x,ix-1)*pow(y,iy)*pow(z,iz-1)*exp0 : 0)
                     + (ix > 0 ? ix*pow(x,ix-1)*pow(y,iy)*pow(z,iz)*expz : 0)
                     + (iz > 0 ? iz*pow(x,ix)*pow(y,iy)*pow(z,iz-1)*expx : 0)
                     + cart*expxz;
        tmp[4][ixyz] = (iy > 0 && iz > 0 ? iy*iz*pow(x,ix)*pow(y,iy-1)*pow(z,iz-1)*exp0 : 0)
                     + (iy > 0 ? iy*pow(x,ix)*pow(y,iy-1)*pow(z,iz)*expz : 0)
                     + (iz > 0 ? iz*pow(x,ix)*pow(y,iy)*pow(z,iz-1)*expy : 0)
                     + cart*expyz;
        tmp[5][ixyz] = (iz > 1 ? iz*(iz-1)*pow(x,ix)*pow(y,iy)*pow(z,iz-2)*exp0 : 0)
                     + (iz > 0 ? 2.0*iz*pow(x,ix)*pow(y,iy)*pow(z,iz-1)*expz : 0)
                     + cart*expzz;
      }
    }
    if (spherical_) {
      carsphlist.carsphfunc_call(index, 1, tmp[0], bxx);
      carsphlist.carsphfunc_call(index, 1, tmp[1], bxy);
      carsphlist.carsphfunc_call(index, 1, tmp[2], byy);
      carsphlist.carsphfunc_call(index, 1, tmp[3], bxz);
      carsphlist.carsphfunc_call(index, 1, tmp[4], byz);
      carsphlist.carsphfunc_call(index, 1, tmp[5], bzz);
    } else {
      copy_n(tmp[0], nxyz, bxx);
      copy_n(tmp[1], nxyz, bxy);
      copy_n(tmp[2], nxyz, byy);
      copy_n(tmp[3], nxyz, bxz);
      copy_n(tmp[4], nxyz, byz);
      copy_n(tmp[5], nxyz, bzz);
    }
    bxx += nxyz;
    bxy += nxyz;
    byy += nxyz;
    bxz += nxyz;
    byz += nxyz;
    bzz += nxyz;
    ++range;
  }
}


shared_ptr<const Shell> Shell::uncontract() const {
  vector<vector<double>> conts;
  vector<pair<int, int>> ranges;
  for (int i = 0; i != exponents_.size(); ++i) {
    vector<double> cont(i, 0.0);
    cont.insert(cont.end(),1.0);
    conts.push_back(cont);
    ranges.push_back({i, i+1});
  }

  auto citer = ranges.begin();
  for (auto iter = conts.begin(); iter != conts.end(); ++iter, ++citer) {
    auto eiter = exponents_.begin();
    double denom = 1.0;
    for (int ii = 2; ii <= angular_number_; ++ii) denom *= 2 * ii - 1;
    for (auto diter = iter->begin(); diter != iter->end(); ++diter, ++eiter)
      *diter *= pow(2.0 * *eiter / pi__, 0.75) * pow(sqrt(4.0 * *eiter), static_cast<double>(angular_number_)) / sqrt(denom);
  }

  return make_shared<Shell>(spherical_, position_, angular_number_, exponents_, conts, ranges);
}
