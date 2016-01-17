//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: petite.cc
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


#include <src/molecule/petite.h>
#include <array>
#include <memory>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace bagel;

static inline array<double,3> matmul33(const vector<double>& a, const array<double,3>& b) {
  assert(a.size() == 9 && b.size() == 3);
  // note that the array is writen in C style (not in fortran style)!!!
  array<double,3> out;
  out[0] = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  out[1] = a[3] * b[0] + a[4] * b[1] + a[5] * b[2];
  out[2] = a[6] * b[0] + a[7] * b[1] + a[8] * b[2];
  return out;
};

Petite::Petite(const vector<shared_ptr<const Atom>>& atoms, const string sym) : sym_(sym) {
  const string c1("c1");
  const string cs("cs");
  const string ci("ci");
  const string c2("c2");
  const string d2("d2");
  const string c2h("c2h");
  const string c2v("c2v");
  const string d2h("d2h");

  natom_ = atoms.size();

  vector<shared_ptr<const Shell>> vbb;
  vector<int> offset;
  int cnt = 0;
  for (auto& a : atoms) {
    vector<shared_ptr<const Shell>> tmp = a->shells();
    vbb.insert(vbb.end(), tmp.begin(), tmp.end());
    offset.push_back(cnt);
    cnt += tmp.size();
  }
  nshell_ = vbb.size();

  if (sym == c1) {
    nirrep_ = 1;
  } else {
    if (sym == c2v){
      SymC2v datc2v;
      symop_ = datc2v.symop();
      nirrep_ = datc2v.nirrep();
    } else if (sym == d2h) {
      SymD2h datd2h;
      symop_ = datd2h.symop();
      nirrep_ = datd2h.nirrep();
    } else if (sym == cs) {
      SymCs datcs;
      symop_ = datcs.symop();
      nirrep_ = datcs.nirrep();
    } else if (sym == ci) {
      SymCi datci;
      symop_ = datci.symop();
      nirrep_ = datci.nirrep();
    } else if (sym == c2) {
      SymC2 datc2;
      symop_ = datc2.symop();
      nirrep_ = datc2.nirrep();
    } else if (sym == d2) {
      SymD2 datd2;
      symop_ = datd2.symop();
      nirrep_ = datd2.nirrep();
    } else if (sym == c2h) {
      SymC2h datc2h;
      symop_ = datc2h.symop();
      nirrep_ = datc2h.nirrep();
    } else {
      assert(false);
    }
    nsymop_ = symop_.size();

    // making map for atoms
    for (int iatom = 0; iatom != natom_; ++iatom) {
      const array<double,3> position = atoms[iatom]->position();
      vector<int> tmp(nsymop_);

      for (int iop = 0; iop != nsymop_; ++iop) {
        const array<double,3> target = matmul33(symop_[iop], position);
        bool found = false;
        for (int jatom = 0; jatom != natom_; ++jatom) {
          const array<double,3> current = atoms[jatom]->position();
          if (current == target) {
            found = true;
            tmp[iop] = jatom;
            break;
          }
        }
        if (!found)
          throw logic_error("Petite constructor error");
      }
      sym_atommap_.push_back(tmp);

      // making map for shells
      for (int i = 0; i != atoms[iatom]->nshell(); ++i) {
        vector<int> stmp(nsymop_);
        for (int iop = 0; iop != nsymop_; ++iop) {
          stmp[iop] = offset[sym_atommap_[iatom][iop]] + i;
        }
        sym_shellmap_.push_back(stmp);
      }
    } // end of atom loop

    // now we determine p1 and p2
    p1_.resize(nshell_);
    lambda_.resize(nshell_ * nshell_);

    fill(p1_.begin(), p1_.end(), 0);

    for (int i = 0; i != nshell_; ++i) {
      bool skipp1 = false;
      for (int iop = 0; iop != nsymop_; ++iop)
        if (sym_shellmap_[i][iop] < i) skipp1 = true;
      if (skipp1) continue;

      p1_[i] = 1;

      for (int j = i; j != nshell_; ++j) {
        const int ij = i * nshell_ + j;
        int nij = 0;

        bool skipp2 = false;
        for (int iop = 0; iop != nsymop_; ++iop) {
          const int ci = sym_shellmap_[i][iop];
          const int cj = sym_shellmap_[j][iop];
          const int cij = ci * nshell_ + cj;
          if (ij == cij) ++nij;
          else if (cij < ij) skipp2 = true;
        }
        if (skipp2) continue;
        lambda_[ij] = (nsymop_ / nij > 0) ? true : false;
      }
    } // end of shell loop

  } // end of if_c1

}


Petite::~Petite() {

}

