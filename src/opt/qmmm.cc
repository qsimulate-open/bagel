//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: qmmm.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <functional>
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/opt/opt.h>

using namespace std;
using namespace bagel;

void QMMM_Tinker::edit_tinker_input(shared_ptr<const Geometry> current) const {
  system("mv -f tinkin.xyz tinkin.xyz.old");
  ifstream fs_tinker_qmmm_old("tinkin.xyz.old");
  ofstream fs_tinker_qmmm("tinkin.xyz");
  string line;
  getline(fs_tinker_qmmm_old, line);
  fs_tinker_qmmm << line << endl;
  {
    int i = 0;
    while (getline(fs_tinker_qmmm_old, line)) {
      stringstream ss(line);
      int dum;
      string atomnm;
      double x,y,z;
      ss >> dum >> atomnm >> x >> y >> z;
      fs_tinker_qmmm << setw(6) << dum << setw(3) << atomnm << setw(20) << setprecision(10) << current->xyz()->element(0, i) * au2angstrom__ <<
        setw(20) << setprecision(10) << current->xyz()->element(1, i) * au2angstrom__ <<
        setw(20) << setprecision(10) << current->xyz()->element(2, i) * au2angstrom__;
      // get the rest dummies and print it to TINKER input
      while (ss >> dum)
        fs_tinker_qmmm << setw(6) << dum;
      fs_tinker_qmmm << endl;
      ++i;
    }
  }
}


void QMMM_Tinker::edit_input(shared_ptr<const Geometry> current) const {
  // generate TINKER input & run TINKER testgrad.x
  // (tinker inputs are prepared in tinker1 and tinker2 directories)

  Timer timer;
  int natom = current->natom();

  chdir("tinker1");

  {
    edit_tinker_input(current);
    system("testgrad -k tinkin.key tinkin.xyz y n n > gradientls");
    system("grep 'Total Potential' gradientls | cut -c 29- > ../energy_qmmm");
    stringstream ss;
    ss << "grep -A " << natom + 3 << " 'Cartesian Gradient Breakdown over Individual Atoms' gradientls | grep -v 'Gradient' | grep 'Anlyt' | cut -c 8- > ../grad_qmmm";
    string pp = ss.str();
    system(pp.c_str());
  }

  chdir("../");
  timer.tick_print("Running TINKER for whole region");
  chdir("tinker2");

  {
    edit_tinker_input(current);
    system("testgrad -k tinkin.key tinkin.xyz y n n > gradientls");
    system("grep 'Total Potential' gradientls | cut -c 29- > ../energy_qm");
    stringstream ss;
    ss << "grep -A " << natom + 3 << " 'Cartesian Gradient Breakdown over Individual Atoms' gradientls | grep -v 'Gradient' | grep 'Anlyt' | cut -c 8- > ../grad_qm";
    string pp = ss.str();
    system(pp.c_str());
  }

  chdir("../");
  timer.tick_print("Running TINKER for QM region");
}


tuple<double,shared_ptr<GradFile>> QMMM_Tinker::do_grad(const int natom) const {
  Timer timer;
  double mmen;
  auto out = make_shared<GradFile>(natom);

  // parse all the outputs
  double mmen_1;
  {
    ifstream fs_energy_qmmm("energy_qmmm");
    if (!fs_energy_qmmm.is_open())
      throw runtime_error("energy_qmmm cannot be opened");
    string line;
    {
      getline(fs_energy_qmmm, line);
      stringstream ss(line);
      ss >> mmen_1;
    }
  }

  double mmen_2;
  {
    ifstream fs_energy_qm("energy_qm");
    if (!fs_energy_qm.is_open())
      throw runtime_error("energy_qm cannot be opened");
    string line;
    {
      getline(fs_energy_qm, line);
      stringstream ss(line);
      ss >> mmen_2;
    }
  }

  auto grad_1 = make_shared<GradFile>(natom);
  {
    ifstream fs_grad_qmmm("grad_qmmm");
    if (!fs_grad_qmmm.is_open())
      throw runtime_error("grad_qmmm cannot be opened");
    string line;
    {
      int i = 0;
      while (getline(fs_grad_qmmm, line)) {
        int dum;
        stringstream ss(line);
        ss >> dum >> grad_1->element(0, i) >> grad_1->element(1, i) >> grad_1->element(2, i);
        ++i;
      }
    }
  }

  auto grad_2 = make_shared<GradFile>(natom);
  {
    ifstream fs_grad_qm("grad_qm");
    if (!fs_grad_qm.is_open())
      throw runtime_error("grad_qm cannot be opened");
    string line;
    {
      int i = 0;
      while (getline(fs_grad_qm, line)) {
        int dum;
        stringstream ss(line);
        ss >> dum >> grad_2->element(0, i) >> grad_2->element(1, i) >> grad_2->element(2, i);
        ++i;
      }
    }
  }

  // energy : kcal / mol
  mmen = (mmen_1 - mmen_2) * kcal2kj__ / (au2kjmol__);

  stringstream ss; ss << "MM energy = " << setw(10) << setprecision(5) << mmen;
  timer.tick_print(ss.str());

  // gradient : kcal / mol / angstrom
  grad_1->scale(kcal2kj__ * au2angstrom__ / au2kjmol__);
  grad_2->scale(kcal2kj__ * au2angstrom__ / au2kjmol__);

  *out = *grad_1 - *grad_2;

  // return the energy and gradient
  return tie(mmen, out);
}
