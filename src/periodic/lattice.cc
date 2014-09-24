//
// BAGEL - Parallel electron correlation program.
// Filename: lattice.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/lattice.h>

using namespace std;
using namespace bagel;

Lattice::Lattice(const shared_ptr<const Geometry> g) : primitive_cell_(g) {

  init();
}

double Lattice::compute_nuclear_repulsion() const {

  double out = 0.0;

  auto cell0 = make_shared<const Geometry>(*primitive_cell_);
  vector<shared_ptr<const Atom>> atoms0 = cell0->atoms();
  int icell0 = 0;
  for (int i = 0; i != ndim_; ++i) icell0 += ncell_ * pow(2 * ncell_ + 1, i);
  int count = 0;
  for (auto& disp : lattice_vectors_) {
    auto cell = make_shared<const Geometry>(*primitive_cell_, disp);
    vector<shared_ptr<const Atom>> atoms = cell->atoms();
    for (auto iter0 = atoms0.begin(); iter0 != atoms0.end(); ++iter0) {
      const double c0 = (*iter0)->atom_charge();
      auto ia0 = distance(atoms0.begin(), iter0);
      for (auto iter1 = atoms.begin(); iter1 != atoms.end(); ++iter1) {
        const double c = (*iter1)->atom_charge();
        auto ia1 = distance(atoms.begin(), iter1);
        if (count == icell0 && ia0 == ia1) continue;
        out += c0 * c / (*iter0)->distance(*iter1);
      }
    }
    ++ count;
  }

  return out;
}

void Lattice::init() {

  ndim_ = primitive_cell_->primitive_vectors().size();
  if (ndim_ > 3)
    cout << "  *** Warning: Dimension in P-SCF is greater than 3!" << endl;
  primitive_rvectors_.resize(ndim_);

  /* TODO: temp parameters */
  ncell_ = 25;
  q_ = 10;

  const int num_vec = pow(2*ncell_+1, ndim_);
  lattice_vectors_.resize(num_vec);
  num_lattice_pts_ = num_vec * primitive_cell_->natom();

  /* Set up lattice vectors */
  switch (ndim_) {
    case 1:
      {
        const array<double, 3> a1 = primitive_cell_->primitive_vectors(0);
        array<double, 3> disp;
        for (int i1 = -ncell_; i1 <= ncell_; ++i1) {
          disp[0] = i1 * a1[0];
          disp[1] = i1 * a1[1];
          disp[2] = i1 * a1[2];
          lattice_vectors_[i1 + ncell_] = disp;
        }
        const double a1sq = dot(a1, a1);
        volume_ = sqrt(a1sq);
        for (int i = 0; i != 3; ++i) primitive_rvectors_[0][i] = 2.0 * pi__ * a1[i] / a1sq;
        break;
      }
    case 2:
      {
        const array<double, 3> a1 = primitive_cell_->primitive_vectors(0);
        const array<double, 3> a2 = primitive_cell_->primitive_vectors(1);
        array<double, 3> disp;
        int count = 0;
        for (int i1 = -ncell_; i1 <= ncell_; ++i1) {
          for (int i2 = -ncell_; i2 <= ncell_; ++i2, ++count) {
            disp[0] = i1 * a1[0] + i2 * a2[0];
            disp[1] = i1 * a1[1] + i2 * a2[1];
            disp[2] = i1 * a1[2] + i2 * a2[2];
            lattice_vectors_[count] = disp;
          }
        }
        array<double, 3> a12 = cross(a1, a2);
        const double a12sq = dot(a12, a12);
        volume_ = sqrt(a12sq);
        const double scale = 2.0 * pi__ / a12sq;
        primitive_rvectors_[0] = cross(a2, a12, scale);
        primitive_rvectors_[1] = cross(a12, a1, scale);
        break;
      }
    case 3:
      {
        const array<double, 3> a1 = primitive_cell_->primitive_vectors(0);
        const array<double, 3> a2 = primitive_cell_->primitive_vectors(1);
        const array<double, 3> a3 = primitive_cell_->primitive_vectors(2);
        array<double, 3> disp;
        int count = 0;
        for (int i1 = -ncell_; i1 <= ncell_; ++i1) {
          for (int i2 = -ncell_; i2 <= ncell_; ++i2) {
            for (int i3 = -ncell_; i3 <= ncell_; ++i3, ++count) {
              disp[0] = i1 * a1[0] + i2 * a2[0] + i3 * a3[0];
              disp[1] = i1 * a1[1] + i2 * a2[1] + i3 * a3[1];
              disp[2] = i1 * a1[2] + i2 * a2[2] + i3 * a3[2];
              lattice_vectors_[count] = disp;
            }
          }
        }
        array<double, 3> a23 = cross(a2, a3);
        volume_ = sqrt(dot(a1, a23));
        const double scale = 2.0 * pi__ / volume_;
        primitive_rvectors_[0] = cross(a2, a3, scale);
        primitive_rvectors_[1] = cross(a3, a1, scale);
        primitive_rvectors_[2] = cross(a1, a2, scale);
        break;
      }
  }

  nuclear_repulsion_ = compute_nuclear_repulsion();
}

double Lattice::dot(array<double, 3> b, array<double, 3> c) { return b[0] * c[0] + b[1] * c[1] + b[2] * c[2]; }

array<double, 3> Lattice::cross(array<double, 3> b, array<double, 3> c, double s) {

  array<double, 3> out;
  out[0] = (b[1] * c[2] - b[2] * c[1]) / s;
  out[1] = (b[2] * c[0] - b[0] * c[2]) / s;
  out[2] = (b[0] * c[1] - b[1] * c[0]) / s;

  return out;
}

void Lattice::generate_kpoints() { /* Monkhorst and Pack PRB 13, 5188 */

  vector<double> u(q_, 0.0);
  int count = 0;
  for (int r = 1; r <= q_; ++r)
    u[count++] = (2.0 * r - q_ - 1.0) / (2.0 * q_);

  const int n = pow(q_, ndim_);
  lattice_rvectors_.resize(n);

  /* set up recriprocal lattice vectors */
  switch(ndim_) {
    case 1:
      {
        const array<double, 3> b1 = primitive_rvectors_[0];
        for (int i = 0; i != q_; ++i) {
          array<double, 3> kvec;
          kvec[0] = u[i] * b1[0];
          kvec[1] = u[i] * b1[1];
          kvec[2] = u[i] * b1[2];
          lattice_rvectors_[i] = kvec;
        }
        break;
      }
    case 2:
      {
        const array<double, 3> b1 = primitive_rvectors_[0];
        const array<double, 3> b2 = primitive_rvectors_[1];
        for (int i = 0; i != q_; ++i) {
          for (int j = 0; j != q_; ++j) {
            array<double, 3> kvec;
            kvec[0] = u[i] * b1[0] + u[j] * b2[0];
            kvec[1] = u[i] * b1[1] + u[j] * b2[1];
            kvec[2] = u[i] * b1[2] + u[j] * b2[2];
            lattice_rvectors_[i * q_ + j] = kvec;
          }
        }
        break;
      }
    case 3:
      {
        const array<double, 3> b1 = primitive_rvectors_[0];
        const array<double, 3> b2 = primitive_rvectors_[1];
        const array<double, 3> b3 = primitive_rvectors_[2];
        for (int i = 0; i != q_; ++i) {
          for (int j = 0; j != q_; ++j) {
            for (int k = 0; k != q_; ++k) {
              array<double, 3> kvec;
              kvec[0] = u[i] * b1[0] + u[j] * b2[0] + u[k] * b3[0];
              kvec[1] = u[i] * b1[1] + u[j] * b2[1] + u[k] * b3[1];
              kvec[2] = u[i] * b1[2] + u[j] * b2[2] + u[k] * b3[2];
              lattice_rvectors_[i * q_ * q_ + j * q_ + k] = kvec;
            }
          }
        }
        break;
      }
  }

}

void Lattice::print_primitive_vectors() const {

  const string indent = "  ";
  cout << indent << "=== Primitive lattice vector(s) ===" << endl << indent << endl;

  for (int i = 0; i != ndim_; ++i)
    cout << indent << fixed << setprecision(6) << "(" << setw(10) << primitive_cell_->primitive_vectors(i)[0] << ", "
                                                      << setw(10) << primitive_cell_->primitive_vectors(i)[1] << ", "
                                                      << setw(10) << primitive_cell_->primitive_vectors(i)[2] << ") " << endl;
  cout << endl;
}

void Lattice::print_lattice_vectors() const {

  const string indent = "  ";
  cout << indent << "=== Lattice vectors ===" << endl << indent << endl;

  for (auto& vec: lattice_vectors_)
    cout << indent << fixed << setprecision(6) << "(" << setw(10) << vec[0] << ", "
                                                      << setw(10) << vec[1] << ", "
                                                      << setw(10) << vec[2] << ") " << endl;
  cout << endl;
}

void Lattice::print_reciprocal_lattice_vectors() const {

  const string indent = "  ";
  cout << indent << "=== Reciprocal Lattice vectors ===" << endl << indent << endl;

  for (auto& vec: lattice_rvectors_)
    cout << indent << fixed << setprecision(6) << "(" << setw(10) << vec[0] << ", "
                                                      << setw(10) << vec[1] << ", "
                                                      << setw(10) << vec[2] << ") " << endl;
  cout << endl;
}

void Lattice::print_lattice_coordinates() const {

  std::ofstream ofs;
  ofs.open("lattice.xyz");
  ofs << num_lattice_pts_ << endl;
  ofs << "[Lattice XYZ Format]" << endl;

  for (auto& disp : lattice_vectors_) {
    auto cell = make_shared<const Geometry>(*primitive_cell_, disp);
    for (auto& atom : cell->atoms()) {
      string name = atom->name();
      name[0] = toupper(name[0]);
      ofs << name << fixed << setprecision(6) << setw(14) << atom->position(0) * au2angstrom__ << "   "
                                              << setw(14) << atom->position(1) * au2angstrom__ << "   "
                                              << setw(14) << atom->position(2) * au2angstrom__ << endl;
    }
  }
}
