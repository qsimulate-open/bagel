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


#include <src/molecule/shell.h>
#include <src/integral/carsph_matrix.h>
#include <src/integral/carsphlist.h>
#include <src/integral/os/overlapbatch.h>
#include <src/integral/os/momentumbatch.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/integral/compos/complexmomentumbatch.h>

using namespace std;
using namespace bagel;

static const CarSphList carsphlist;

Shell::Shell(const bool sph, const array<double,3>& _position, int _ang, const vector<double>& _expo,
                       const vector<vector<double>>& _contr,  const vector<pair<int, int>>& _range, const array<double,3>& _vector_potential)
 : Shell_base(sph, _position, _ang),
   exponents_(_expo), contractions_(_contr), contraction_ranges_(_range), dummy_(false), relativistic_(false), vector_potential_(_vector_potential) {

#if 0
  vector_potential_ = {{ 0.0, 0.0, 0.0 }};
  cout << "Caution:  Using Gaussian orbitals with common gauge origin." << endl;
#endif
#if 0
  vector_potential_ = {{ 0.072, -0.05, 0.18 }};
  cout << "Caution:  Using a constant phase factor in Shell - for debugging only." << endl;
  // Remember to also turn on the version for the dummy atoms.
#endif

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


Shell::Shell(const bool sph) : Shell_base(sph), exponents_{0.0}, contractions_{{1.0}},
                               contraction_ranges_{make_pair(0,1)}, dummy_(true), vector_potential_{{0.0,0.0,0.0}} {
#if 0
  vector_potential_ = {{ 0.072, -0.05, 0.18 }};
  cout << "Caution:  Using a constant phase factor in Shell - for debugging only." << endl;
#endif
  contraction_lower_.push_back(0);
  contraction_upper_.push_back(1);
}


std::shared_ptr<const Shell> Shell::move_atom(const array<double,3>& displacement) const {
  auto out = make_shared<Shell>(*this);
  out->position_[0] += displacement[0];
  out->position_[1] += displacement[1];
  out->position_[2] += displacement[2];
  return out;
}


std::shared_ptr<const Shell> Shell::move_atom(const double* displacement) const {
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
        range.push_back(make_pair(contraction_ranges_[i].first-smallest, contraction_ranges_[i].second-smallest));
      }
      out.push_back(make_shared<const Shell>(spherical_, position_, angular_number_, expo, contr, range, vector_potential_));
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
  static_assert(increment==1||increment==-1, "illegal call of Shell::kinetic_balance_uncont");
  int i = 0;
  vector<vector<double>> conts;
  vector<pair<int, int>> ranges;
  for (auto e = exponents_.begin(); e != exponents_.end(); ++e, ++i) {
    vector<double> cont(exponents_.size(), 0);
    cont[i] = 1.0;
    conts.push_back(cont);
    ranges.push_back(make_pair(i,i+1));
  }
  return angular_number_+increment < 0 ? nullptr : make_shared<const Shell>(false, position_, angular_number_+increment, exponents_, conts, ranges, vector_potential_);
}

shared_ptr<const Shell> Shell::cartesian_shell() const {
  auto out = make_shared<Shell>(false, position_, angular_number_, exponents_, contractions_, contraction_ranges_, vector_potential_);
  return out;
}

void Shell::init_relativistic() {
  if (angular_number_ == 6) throw runtime_error("Relativistic codes cannot use i-type main basis functions, since j-type would be needed for the small component.");
  relativistic_ = true;
  aux_decrement_ = kinetic_balance_uncont<-1>();
  aux_increment_ = kinetic_balance_uncont<1>();

  // small is a transformation matrix (x,y,z components)
#if 0
  // overlap = S^-1 between auxiliary functions
  shared_ptr<const Matrix> overlap = overlap_compute_();
  small_ = moment_compute(overlap);
  //array<shared_ptr<const Matrix>,3> comparison = new_moment_compute();
#else
  small_ = new_moment_compute();
  //array<shared_ptr<const Matrix>,3> comparison = moment_compute(overlap);
#endif

  /*
  array<shared_ptr<const Matrix>,3> difference;
  for (int i=0; i!=3; i++) difference[i] = make_shared<Matrix>( *(comparison[i]) - *(small_[i]) );

  cout << scientific << "angular_number = " << angular_number_ << ", nprim = " << num_primitive() << ", ncont = " << num_contracted() << endl;
  cout << "x error = " << difference[0]->rms() << ", relative error = " << (difference[0]->rms()/small_[0]->rms()) << endl;
  cout << "y error = " << difference[1]->rms() << ", relative error = " << (difference[1]->rms()/small_[1]->rms()) << endl;
  cout << "z error = " << difference[2]->rms() << ", relative error = " << (difference[2]->rms()/small_[2]->rms()) << endl;
  */
}


void Shell::init_relativistic_london(const array<double,3> magnetic_field) {
  if (angular_number_ == 6) throw runtime_error("Relativistic codes cannot use i-type main basis functions, since j-type would be needed for the small component.");
  relativistic_ = true;
  aux_decrement_ = kinetic_balance_uncont<-1>();
  aux_increment_ = kinetic_balance_uncont<1>();

  // small is a transformation matrix (x,y,z components)
#if 0
  // overlap = S^-1 between auxiliary functions
  shared_ptr<const ZMatrix> overlap = make_shared<const ZMatrix>(*overlap_compute_(), 1.0);
  zsmall_ = moment_compute(overlap, magnetic_field);
  for (int i=0; i!=3; i++) zsmallc_[i] = zsmall_[i]->get_conjg();
  //array<shared_ptr<const ZMatrix>,3> comparison = new_moment_compute(magnetic_field);
#else
  zsmall_ = new_moment_compute(magnetic_field);
  for (int i=0; i!=3; i++) zsmallc_[i] = zsmall_[i]->get_conjg();
  //array<shared_ptr<const ZMatrix>,3> comparison = moment_compute(overlap, magnetic_field);
#endif

  /*
  array<shared_ptr<const ZMatrix>,3> difference;
  for (int i=0; i!=3; i++) difference[i] = make_shared<ZMatrix>( *(comparison[i]) - *(zsmall_[i]) );

  cout << scientific << "angular_number = " << angular_number_ << ", nprim = " << num_primitive() << ", ncont = " << num_contracted() << endl;
  cout << "x error = " << difference[0]->rms() << ", relative error = " << (difference[0]->rms()/zsmall_[0]->rms()) << endl;
  cout << "y error = " << difference[1]->rms() << ", relative error = " << (difference[1]->rms()/zsmall_[1]->rms()) << endl;
  cout << "z error = " << difference[2]->rms() << ", relative error = " << (difference[2]->rms()/zsmall_[2]->rms()) << endl;
  */

  // to turn off magnetic balance
  //array<shared_ptr<const Matrix>,3> moment = moment_compute(overlap_compute_());
  //for (int i=0; i!=3; i++) { zsmall_[i] = make_shared<ZMatrix>(*moment[i], 1.0); }
  //for (int i=0; i!=3; i++) zsmallc_[i] = zsmall_[i]->get_conjg();

}


shared_ptr<const Matrix> Shell::overlap_compute_() const {

  const int asize_inc = aux_increment_->nbasis();
  const int asize_dec = aux_decrement_ ? aux_decrement_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  auto overlap = make_shared<Matrix>(a,a, true);

  {
    OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_increment_, aux_increment_}});
    ovl.compute();
    for (int i = 0; i != aux_increment_->nbasis(); ++i)
      copy_n(ovl.data() + i*asize_inc, asize_inc, overlap->element_ptr(0,i));
  }
  if (aux_decrement_) {
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_decrement_, aux_decrement_}});
      ovl.compute();
      for (int i = 0; i != asize_dec; ++i)
        copy_n(ovl.data() + i*asize_dec, asize_dec, overlap->element_ptr(asize_inc, i+asize_inc));
    }
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_increment_, aux_decrement_}});
      ovl.compute();
      for (int i = 0; i != asize_dec; ++i)
        for (int j = 0; j != asize_inc; ++j)
          overlap->element(j,i+asize_inc) = overlap->element(i+asize_inc,j) = *(ovl.data()+j+asize_inc*i);
    }
  }

  return overlap;
}


array<shared_ptr<const Matrix>,3> Shell::moment_compute(const shared_ptr<const Matrix> overlap) const {
  const int ssize = nbasis();
  const int asize_inc = aux_increment_->nbasis();
  const int asize_dec = aux_decrement_ ? aux_decrement_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  auto coeff0 = make_shared<MomentumBatch>(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_increment_}});
  coeff0->compute();

  shared_ptr<MomentumBatch> coeff1;
  if (aux_decrement_) {
    coeff1 = make_shared<MomentumBatch>(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_decrement_}});
    coeff1->compute();
  } else {
    // just to run. coeff1 is not referenced in the code
    coeff1 = coeff0;
  }

  const double* carea0 = coeff0->data();
  const double* carea1 = coeff1->data();

  auto tmparea = make_shared<Matrix>(ssize,a, true);
  array<shared_ptr<const Matrix>,3> out;

  const static CarSphList carsphlist;
  for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
    if (spherical_) {
      const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
      const int nloop = num_contracted() * asize_inc;
      carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea->data());
    } else {
      const size_t size = coeff0->asize_final()*num_contracted()*num_primitive();
      assert(size == asize_inc*ssize);
      copy(carea0, carea0+size, tmparea->data());
    }
    if (aux_decrement_) {
      if (spherical_) {
        const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
        const int nloop = num_contracted() * asize_dec;
        carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea->data()+asize_inc*ssize);
      } else {
        const size_t size = coeff1->asize_final()*num_contracted()*num_primitive();
        assert(size == asize_dec*ssize);
        copy(carea1, carea1+size, tmparea->data()+asize_inc*ssize);
        }
    }

    out[i] = tmparea->transpose()->solve(overlap, overlap->ndim());
  }
  return out;
}


array<shared_ptr<const ZMatrix>,3> Shell::moment_compute(const shared_ptr<const ZMatrix> overlap, const array<double, 3> magnetic_field) const {
  const int ssize = nbasis();
  const int asize_inc = aux_increment_->nbasis();
  const int asize_dec = aux_decrement_ ? aux_decrement_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  auto coeff0 = make_shared<ComplexMomentumBatch>(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_increment_}}, magnetic_field);
  coeff0->compute();

  shared_ptr<ComplexMomentumBatch> coeff1;
  if (aux_decrement_) {
    coeff1 = make_shared<ComplexMomentumBatch>(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_decrement_}}, magnetic_field);
    coeff1->compute();
  } else {
    // just to run. coeff1 is not referenced in the code
    coeff1 = coeff0;
  }

  const complex<double>* carea0 = coeff0->data();
  const complex<double>* carea1 = coeff1->data();

  auto tmparea = make_shared<ZMatrix>(ssize,a, true);
  array<shared_ptr<const ZMatrix>,3> out;

  const static CCarSphList carsphlist;
  for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
    if (spherical_) {
      const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
      const int nloop = num_contracted() * asize_inc;
      carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea->data());
    } else {
      const size_t size = coeff0->asize_final()*num_contracted()*num_primitive();
      assert(size == asize_inc*ssize);
      copy(carea0, carea0+size, tmparea->data());
    }
    if (aux_decrement_) {
      if (spherical_) {
        const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
        const int nloop = num_contracted() * asize_dec;
        carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea->data()+asize_inc*ssize);
      } else {
        const size_t size = coeff1->asize_final()*num_contracted()*num_primitive();
        assert(size == asize_dec*ssize);
        copy(carea1, carea1+size, tmparea->data()+asize_inc*ssize);
        }
    }

    const complex<double> ii (0.0, 1.0);
    tmparea->scale(ii);

    out[i] = tmparea->transpose_conjg()->solve(overlap, overlap->ndim());

  }
  return out;
}


// In DFT we want to compute values of basis functions on grid
void Shell::compute_grid_value(double* b, double* dx, double* dy, double* dz, const double& x, const double& y, const double& z) const {
  const double rr = x*x+y*y+z*z;
  auto range = contraction_ranges_.begin();
  const int nang = angular_number();

  double tmp0[65];
  double tmpx[65];
  double tmpy[65];
  double tmpz[65];
  static_assert(65 > (ANG_HRR_END+1)*(ANG_HRR_END+1), "ANG_HRR_END is assumed to be 7");
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
    if (fabs(exp0) >= 1.0e-14) {
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
    }
    b  += nxyz;
    dx += nxyz;
    dy += nxyz;
    dz += nxyz;
    ++range;
  }

}


// In DFT we want to compute values of basis functions on grid
void Shell::compute_grid_value_deriv2(double* bxx, double* bxy, double* byy, double* bxz, double* byz, double* bzz,
                                      const double& x, const double& y, const double& z) const {
  const double rr = x*x+y*y+z*z;
  auto range = contraction_ranges_.begin();
  double tmp[6][50];
  assert(50 > ANG_HRR_END*ANG_HRR_END);

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


// for each primitive basis function
array<shared_ptr<const Matrix>,6> Shell::mblock(const double exponent) const {

  const int norig = (angular_number_+1) * (angular_number_+2) / 2;
  const int ninc = norig + angular_number_ + 2;
  const int ndec = norig - angular_number_ - 1;
  assert(ninc == (aux_increment_->angular_number_+1) * (aux_increment_->angular_number_+2) / 2);
  if (aux_decrement_) assert(ndec == (aux_decrement_->angular_number_+1) * (aux_decrement_->angular_number_+2) / 2);
  else assert (ndec == 0);

  array<shared_ptr<Matrix>,6> mcart;
  for (int i=0; i!=3; i++) mcart[i] = make_shared<Matrix>(ninc, norig, true);
  if (aux_decrement_) for (int i=3; i!=6; i++) mcart[i] = make_shared<Matrix>(ndec, norig, true);
  for (int i=0; i!=norig; i++) {
    int x, y, z;
    size_t column = 0;
    for (int i=0; i<=angular_number_; i++) {
      z = i;
      for (int j=0; j<=(angular_number_-i); j++) {
        y = j;
        x = angular_number_ - i - j;

        // three components of the angular momentum
        array<int,3> index = {{x, y, z}};
        array<double,3> dindex;
        for (int k=0; k!=3; k++) dindex[k] = index[k];

        assert(column == index[2]*(angular_number_+1) - index[2]*(index[2]-1)/2 + index[1]);
        const double talph = 2.0 * exponent;

        // k tells us which dimension of the momentum operator we're using
        for (int k=0; k!=3; k++) {

          // -i a_x phi^(x-1)
          if (index[k]!=0) {
            array<int,3> newindex = index;
            newindex[k]--;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ - 1);
            const size_t row = newindex[2]*(angular_number_) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[3+k]->element(row, column) = -dindex[k];
          }

          // +i 2alpha phi^(x+1)
          {
            array<int,3> newindex = index;
            newindex[k]++;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ + 1);
            const size_t row = newindex[2]*(angular_number_+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[k]->element(row, column) = talph;
          }
        }
        column++;
      }
    }
  }

  array<shared_ptr<const Matrix>,6> out;
  const int nblock = aux_decrement_ ? 6 : 3;

  // convert this block from cartesian to spherical
  if (spherical_) for (int i=0; i!=nblock; i++) out[i] = make_shared<const Matrix>(*mcart[i] * *carsph_matrix(angular_number_));
  else for (int i=0; i!=nblock; i++) out[i] = mcart[i];

  return out;
}


// for each primitive basis function
array<shared_ptr<const ZMatrix>,6> Shell::mblock(const double exponent, const array<double,3> magnetic_field) const {

  const int norig = (angular_number_+1) * (angular_number_+2) / 2;
  const int ninc = norig + angular_number_ + 2;
  const int ndec = norig - angular_number_ - 1;
  const complex<double> imag (0.0, 1.0);
  assert(ninc == (aux_increment_->angular_number_+1) * (aux_increment_->angular_number_+2) / 2);
  if (aux_decrement_) assert(ndec == (aux_decrement_->angular_number_+1) * (aux_decrement_->angular_number_+2) / 2);
  else assert (ndec == 0);

  array<shared_ptr<ZMatrix>,6> mcart;
  for (int i=0; i!=3; i++) mcart[i] = make_shared<ZMatrix>(ninc, norig, true);
  if (aux_decrement_) for (int i=3; i!=6; i++) mcart[i] = make_shared<ZMatrix>(ndec, norig, true);
  for (int i=0; i!=norig; i++) {
    int x, y, z;
    size_t column = 0;
    for (int i=0; i<=angular_number_; i++) {
      z = i;
      for (int j=0; j<=(angular_number_-i); j++) {
        y = j;
        x = angular_number_ - i - j;

        // three components of the angular momentum
        array<int,3> index = {{x, y, z}};
        array<double,3> dindex;
        for (int k=0; k!=3; k++) dindex[k] = index[k];
        const array<const int,3> fwd  = {{1, 2, 0}};
        const array<const int,3> back = {{2, 0, 1}};

        assert(column == index[2]*(angular_number_+1) - index[2]*(index[2]-1)/2 + index[1]);
        const complex<double> tialph = imag * 2.0 * exponent;
        const array<const complex<double>,3> halfb = {{0.5*magnetic_field[0], 0.5*magnetic_field[1], 0.5*magnetic_field[2]}};

        // k tells us which dimension of the momentum operator we're using
        for (int k=0; k!=3; k++) {

          // -i a_x phi^(x-1)
          if (index[k]!=0) {
            array<int,3> newindex = index;
            newindex[k]--;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ - 1);
            const size_t row = newindex[2]*(angular_number_) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[3+k]->element(row, column) = -dindex[k]*imag;
          }

          // +i 2alpha phi^(x+1)
          {
            array<int,3> newindex = index;
            newindex[k]++;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ + 1);
            const size_t row = newindex[2]*(angular_number_+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[k]->element(row, column) = tialph;
          }

          // + 1/2 B_y phi^(z+1)
          {
            array<int,3> newindex = index;
            newindex[back[k]]++;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ + 1);
            const size_t row = newindex[2]*(angular_number_+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[k]->element(row, column) = halfb[fwd[k]];
          }

          // - 1/2 B_z phi^(y+1)
          {
            array<int,3> newindex = index;
            newindex[fwd[k]]++;
            assert(newindex[0] + newindex[1] + newindex[2] == angular_number_ + 1);
            const size_t row = newindex[2]*(angular_number_+2) - newindex[2]*(newindex[2]-1)/2 + newindex[1];
            mcart[k]->element(row, column) = -halfb[back[k]];
          }
        }
        column++;
      }
    }
  }

  array<shared_ptr<const ZMatrix>,6> out;
  const int nblock = aux_decrement_ ? 6 : 3;
  for (int i=0; i!=nblock; i++) mcart[i]->scale(complex<double>(0.0, -1.0));

  // convert this block from cartesian to spherical
  if (spherical_) for (int i=0; i!=nblock; i++) out[i] = make_shared<const ZMatrix>(*mcart[i] * *make_shared<ZMatrix>(*carsph_matrix(angular_number_), 1.0));
  else for (int i=0; i!=nblock; i++) out[i] = mcart[i];

  return out;
}


array<shared_ptr<const Matrix>,3> Shell::new_moment_compute() const {

  int norig;
  int n;
  norig = ( (angular_number_+1) * (angular_number_+2) / 2 );
  const int ninc = norig + (angular_number_ + 2);
  const int ndec = norig - (angular_number_ + 1);
  n = ninc + ndec;
  if (spherical_ ) norig = 2*angular_number_+1;

  assert(ninc == (aux_increment_->angular_number_+1) * (aux_increment_->angular_number_+2) / 2 );
  if (aux_decrement_) assert(ndec == (aux_decrement_->angular_number_+1) * (aux_decrement_->angular_number_+2) / 2 );
  else assert (ndec == 0);
  assert(aux_increment_->num_primitive() == num_primitive());
  if(aux_decrement_)assert(aux_decrement_->num_primitive() == num_primitive());

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2 blocks, one for L+1 & one for L-1
  array<shared_ptr<Matrix>,3> tmp;
  for (int i=0; i!=3; i++) tmp[i] = make_shared<Matrix>(n*num_primitive(), norig*num_primitive(), true);
  for (int j=0; j!=num_primitive(); j++) {
    auto thisblock = mblock(exponents_[j]);
    for (int i=0; i!=3; i++) tmp[i]->copy_block( j*ninc, j*norig, ninc, norig, thisblock[i] );
    if (aux_decrement_) for (int i=0; i!=3; i++) tmp[i]->copy_block( ninc*num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i] );
  }

  // build the contraction matrix
  auto contract = make_shared<Matrix>(norig*num_primitive(), norig*num_contracted(), true);
  for (int j=0; j!=num_contracted(); j++) {
    for (int k=contraction_ranges_[j].first; k!=contraction_ranges_[j].second; k++) {
      for (int l=0; l!=norig; l++) {
        contract->element(k*norig+l, j*norig+l) = (contractions_[j])[k];
      }
    }
  }

  // contract
  array<shared_ptr<const Matrix>,3> out;
  for (int i=0; i!=3; i++) out[i] = make_shared<const Matrix>(*tmp[i] * *contract);
  return out;
}


array<shared_ptr<const ZMatrix>,3> Shell::new_moment_compute(const array<double,3> magnetic_field) const {

  int norig;
  int n;
  norig = ( (angular_number_+1) * (angular_number_+2) / 2 );
  const int ninc = norig + (angular_number_ + 2);
  const int ndec = norig - (angular_number_ + 1);
  if (spherical_ ) norig = 2*angular_number_+1;
  n = ninc + ndec;

  assert(ninc == (aux_increment_->angular_number_+1) * (aux_increment_->angular_number_+2) / 2 );
  if (aux_decrement_) assert(ndec == (aux_decrement_->angular_number_+1) * (aux_decrement_->angular_number_+2) / 2 );
  else assert (ndec == 0);
  assert(aux_increment_->num_primitive() == num_primitive());
  if(aux_decrement_)assert(aux_decrement_->num_primitive() == num_primitive());

  // build the momentum transformation matrix for primitive functions
  // each exponent gets 2 blocks, one for L+1 & one for L-1
  array<shared_ptr<ZMatrix>,3> tmp;
  for (int i=0; i!=3; i++) tmp[i] = make_shared<ZMatrix>(n*num_primitive(), norig*num_primitive(), true);
  for (int j=0; j!=num_primitive(); j++) {
    auto thisblock = mblock(exponents_[j], magnetic_field);
    for (int i=0; i!=3; i++) tmp[i]->copy_block( j*ninc, j*norig, ninc, norig, thisblock[i] );
    if (aux_decrement_) for (int i=0; i!=3; i++) tmp[i]->copy_block( ninc*num_primitive()+j*ndec, j*norig, ndec, norig, thisblock[3+i] );
  }

  // build the contraction matrix
  auto contract = make_shared<ZMatrix>(norig*num_primitive(), norig*num_contracted(), true);
  for (int j=0; j!=num_contracted(); j++) {
    for (int k=contraction_ranges_[j].first; k!=contraction_ranges_[j].second; k++) {
      for (int l=0; l!=norig; l++) {
        contract->element(k*norig+l, j*norig+l) = (contractions_[j])[k];
      }
    }
  }

  // contract
  array<shared_ptr<const ZMatrix>,3> out;
  for (int i=0; i!=3; i++) out[i] = make_shared<const ZMatrix>(*tmp[i] * *contract);
  return out;
}
