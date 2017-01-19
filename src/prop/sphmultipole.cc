//
// BAGEL - Parallel electron correlation program.
// Filename: sphmultipole.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/prop/sphmultipole.h>

using namespace std;
using namespace bagel;

SphMultipole::SphMultipole(shared_ptr<const Geometry> g, shared_ptr<const Matrix> d, const int rank)
  : geom_(g), density_(d), rank_(rank) {

  if (rank_ > 3)
    throw logic_error("Higher-order multipole moments not available");
}

vector<complex<double>> SphMultipole::compute() const { // slow

  const size_t size = (rank_+1)*(rank_+1);
  vector<complex<double>> out(size);

  auto o0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    auto o1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++o1) {

      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          MultipoleBatch mpole(input, geom_->charge_center(), rank_);
          mpole.compute();
          assert(mpole.num_blocks() == size);

          const int dimb1 = input[0]->nbasis();
          const int dimb0 = input[1]->nbasis();
          vector<const complex<double>*> dat(mpole.num_blocks());
          for (int i = 0; i != mpole.num_blocks(); ++i)
            dat[i] = mpole.data() + mpole.size_block()*i;

          for (int i = *offset0; i != dimb0 + *offset0; ++i) {
            for (int j = *offset1; j != dimb1 + *offset1; ++j) {
              for (int k = 0; k != mpole.num_blocks(); ++k) {
                out[k] += *dat[k]++ * density_->element(j,i);
              }
            }
          }

        }
      }
    }
  }

  // Rlmc = (-1)^m * Rlm + Rl-m     iRlms = (-1)^m * Rlm - Rl-m        Rlm = out[l*l+l+m]
  cout << "    * Dipole moment:" << endl;
  const double mpole11c = (-1.0 * out[3] + out[1]).real();
  const double mpole11s = (-1.0 * out[3] - out[1]).imag();
  const double mpole10 = out[2].real();
  cout << "           (" << setw(12) << setprecision(6) << mpole11c << ", " << setw(12)
                                                        << mpole11s << ", " << setw(12)
                                                        << mpole10 << ") a.u." << endl << endl;

  if (rank_ > 1) {
    // nuclear contribution
    array<double,5> qm;
    fill(qm.begin(), qm.end(), 0.0);
    for (auto& a : geom_->atoms()) {
      array<double, 3> rc;
      double rr = 0.0;
      for (int i = 0; i != 3; ++i) {
        rc[i] = a->position(i) - geom_->charge_center()[i];
        rr += rc[i]*rc[i];
      }
      qm[0] += a->atom_charge() * (1.5*rc[0]*rc[0] - 0.5*rr); //xx
      qm[1] += a->atom_charge() * (1.5*rc[1]*rc[1] - 0.5*rr); //yy
      qm[2] += a->atom_charge() * 1.5*rc[0]*rc[1];            //xy
      qm[3] += a->atom_charge() * 1.5*rc[0]*rc[2];            //xz
      qm[4] += a->atom_charge() * 1.5*rc[1]*rc[2];            //yz
    }

    // Stone's convention
    const double sqrt3 = sqrt(3.0);
    const double sqrt3inv = 1.0/sqrt3;
    const double mpole22c = sqrt3inv*(qm[0]-qm[1]) - 2.0*sqrt3*(out[8] + out[4]).real();
    const double mpole22s = 2.0*sqrt3inv*qm[2] - 2.0*sqrt3*(out[8] - out[4]).imag();
    const double mpole21c = 2.0*sqrt3inv*qm[3] - sqrt3*(-1.0 * out[7] + out[5]).real();
    const double mpole21s = 2.0*sqrt3inv*qm[4] - sqrt3*(-1.0 * out[7] - out[5]).imag();
    const double mpole20  = -1.0*(qm[0]+qm[1]) - 2.0*out[6].real();
    cout << "    * Quadrupole moment:" << endl;
    cout << "           (" << setw(12) << setprecision(6);
    if (abs(mpole20)  > 1e-7) cout << "  Q20 = " << setw(11) << mpole20  << ",";
    if (abs(mpole21c) > 1e-7) cout << " Q21c = " << setw(11) << mpole21c << ",";
    if (abs(mpole21s) > 1e-7) cout << " Q21s = " << setw(11) << mpole21s << ",";
    if (abs(mpole22c) > 1e-7) cout << " Q22c = " << setw(11) << mpole22c << ",";
    if (abs(mpole22s) > 1e-7) cout << " Q22s = " << setw(11) << mpole22s;
    cout << ") a.u." << endl;
  }

  if (rank_ > 2) {
    // nuclear contribution
    array<double,7> qm;
    fill(qm.begin(), qm.end(), 0.0);
    for (auto& a : geom_->atoms()) {
      array<double, 3> rc;
      double rr = 0.0;
      for (int i = 0; i != 3; ++i) {
        rc[i] = a->position(i) - geom_->charge_center()[i];
        rr += rc[i]*rc[i];
      }
      qm[0] += a->atom_charge() * (2.5*rc[0]*rc[0]*rc[0] - 1.5*rr*rc[0]); //xxx
      qm[2] += a->atom_charge() * (2.5*rc[0]*rc[1]*rc[1] - 0.5*rr*rc[0]); //xyy
      qm[1] += a->atom_charge() * (2.5*rc[0]*rc[0]*rc[1] - 0.5*rr*rc[1]); //yxx
      qm[3] += a->atom_charge() * (2.5*rc[1]*rc[1]*rc[1] - 1.5*rr*rc[1]); //yyy
      qm[4] += a->atom_charge() * (2.5*rc[0]*rc[0]*rc[2] - 0.5*rr*rc[2]); //zxx
      qm[6] += a->atom_charge() * (2.5*rc[1]*rc[1]*rc[2] - 0.5*rr*rc[2]); //zyy
      qm[5] += a->atom_charge() *  2.5*rc[0]*rc[1]*rc[2];                 //xyz
    }

    // Stone's convention
    const double mpole33c = sqrt(0.1)*(qm[0]-3.0*qm[2]) - 6.0*sqrt(10.0)*(-1.0*out[15]+out[9]).real();
    const double mpole33s = sqrt(0.1)*(3.0*qm[1]-qm[3]) - 6.0*sqrt(10.0)*(-1.0*out[15]-out[9]).imag();
    const double mpole32c = sqrt(0.6)*(qm[4]-qm[6]) - 2.0*sqrt(15.0)*(out[14]+out[10]).real();
    const double mpole32s = 2.0*sqrt(0.6)*qm[5] - 2.0*sqrt(15.0)*(out[14]-out[10]).imag();
    const double mpole31c = -sqrt(1.5)*(qm[0]+qm[2]) - 2.0*sqrt(6.0)*(-1.0*out[13]+out[11]).real();
    const double mpole31s = -sqrt(1.5)*(qm[1]+qm[3]) - 2.0*sqrt(6.0)*(-1.0*out[13]-out[11]).imag();
    const double mpole30  = -1.0*(qm[4]+qm[6]) - 6.0*out[12].real();
    cout << "    * Octopole moment:" << endl;
    cout << "           (" << setw(12) << setprecision(6);
    if (abs(mpole30)  > 1e-7) cout << "  Q30 = " << setw(11) << mpole30  << ",";
    if (abs(mpole31c) > 1e-7) cout << " Q31c = " << setw(11) << mpole31c << ",";
    if (abs(mpole31s) > 1e-7) cout << " Q31s = " << setw(11) << mpole31s << ",";
    if (abs(mpole32c) > 1e-7) cout << " Q32c = " << setw(11) << mpole32c << ",";
    if (abs(mpole32s) > 1e-7) cout << " Q32s = " << setw(11) << mpole32s << ",";
    if (abs(mpole33c) > 1e-7) cout << " Q33c = " << setw(11) << mpole33c << ",";
    if (abs(mpole33s) > 1e-7) cout << " Q33s = " << setw(11) << mpole33s;
    cout << ") a.u." << endl;
  }

  cout << "      about the centre of charge:      (" << geom_->charge_center()[0] << ", " << setw(12)
                                                     << geom_->charge_center()[1] << ", " << setw(12)
                                                     << geom_->charge_center()[2] << ")" << endl;

  return out;
}
