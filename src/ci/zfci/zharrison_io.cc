//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison_io.cc
// Copyright (C) 2016 Toru Shiozaki
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

#include <src/ci/zfci/zharrison.h>

using namespace std;
using namespace bagel;

// The following code is due to George Booth

// Write relativistic integrals to an file in chemists notation suitable for reading into NECI
void ZHarrison::dump_ints() const {

  cout << "" << endl;
  cout << "Writing relativistic integrals to FCIDUMP.REL file" << endl;
  cout << setprecision(10) << endl;
  // kramers_coeff(0) accesses the + functions
  cout << "Number of AO basis functions: " << jop_->coeff()->ndim() << endl; // ndim is length of first index
  cout << "Number of MO basis functions to dump: " << jop_->coeff()->mdim()/2 << endl; //mdim for length for second index

  ofstream fs("FCIDUMP");

  int jfac;
  int j2fac;
  int kfac;
  int k2fac;
  complex<double> val;
  complex<double> tval;

  if (fs.is_open()) {

    fs << " &FCI NORB= " << norb_*2 << ",NELEC= " << nele_ << ",ORBSYM= ";
    for (int i = 0; i != norb_*2; ++i)
      fs << "1, ";
    fs << endl;
    fs << " ISYM= 0 ," << endl;
    fs << " TREL=.TRUE." << endl;
    fs << " &END" << endl;

    fs << setw(20) << setprecision(15);

    for (int i = 0; i != 16; ++i) {
      if (!jop_->mo2e()->exist(i)) continue;
      cout << "Writing 2e integral block " << i+1 << " / 16 : ";
      shared_ptr<const ZMatrix> tmp = jop_->mo2e(i);
      // assuming here that the fastest bit in i corresponds to the slowest orbital
      // in mo2e
      if ((i & 1) == 1) { jfac = 1; } // j is a plus kramers spinor
      else jfac = 2;              // j is a minus kramers spinor
      if ((i & 2) == 2) { j2fac = 1; }
      else j2fac = 2;
      if ((i & 4) == 4) { kfac = 1; }
      else kfac = 2;
      if((i & 8) == 8) { k2fac = 1; }
      else k2fac = 2;

      if ( k2fac == 1) cout << "( + ";
      else cout << "( - ";
      if ( j2fac == 1) cout << "+ |";
      else cout << "- |";
      if ( kfac == 1) cout << " + ";
      else cout << " - ";
      if ( jfac == 1) cout << "+ )" << endl;
      else cout << "- )" << endl;

      for (int j = 0; j != norb_; ++j)
        for (int j2 = 0; j2 != norb_; ++j2)
          for (int k = 0; k != norb_; ++k)
            for (int k2 = 0; k2 != norb_; ++k2) {
              val = tmp->element(k2+norb_*k, j2+norb_*j);
              if (abs(val) > 1.0e-9) {
                fs << setw(20) << val << setw(4) << 2*k2+k2fac << setw(4)
                  << 2*j2+j2fac << setw(4) << 2*k+kfac << setw(4) << 2*j+jfac << endl;   // electron 1, electron 2
              }
            }
    }

    tval = 0.0;
    for (int i = 0; i != 4; ++i) {
      if (!jop_->mo1e()->exist(i)) continue;
      cout << "Writing 1e integral block " << i+1 << " / 4" << endl;
      shared_ptr<const ZMatrix> tmp = jop_->mo1e(i);

      if ((i & 1) == 1) {jfac = 1; }
      else jfac = 2;
      if ((i & 2) == 2) {kfac = 1; }
      else kfac = 2;

      for (int j = 0; j != norb_; ++j) {
        for (int k = 0; k != norb_; ++k) {
          val = tmp->element(k, j);
          if (abs(val) > 1.0e-9)
            fs << val << setw(4) << 2*k+kfac << setw(4) << 2*j+jfac << "   0   0" << endl;
          if ((j == k) && ((i == 0) || (i == 3))) tval += val;
        }
      }
    }
    fs << "(" << jop_->core_energy()  + geom_->nuclear_repulsion() << ",0.0)" << "   0   0   0   0" << endl;
    fs.close();
  }
  else throw runtime_error("Unable to open file");
}


void ZHarrison::read_external_rdm12_av(const string& file) {
  rdm1_av_expanded_ = make_shared<ZRDM<1>>(norb_*2);
  rdm2_av_expanded_ = make_shared<ZRDM<2>>(norb_*2);

  // feed RDM1
  {
    ifstream fs(file + ".rdm1"); 
    if (!fs.is_open()) throw runtime_error(file + ".rdm1 cannot be opened");
    string line;
    while (getline(fs, line)) {
      stringstream ss(line);
      int i, j;
      double re, im;
      ss >> i >> j >> re >> im;
      const int ii = ((i-1)/2) + ((i-1)%2)*norb_;
      const int jj = ((j-1)/2) + ((j-1)%2)*norb_;
      const complex<double> dat(re, im);
      rdm1_av_expanded_->element(ii, jj) = dat; 
      rdm1_av_expanded_->element(jj, ii) = conj(dat);
    }
    fs.close();
  }
  {
    ifstream fs(file + ".rdm2"); 
    if (!fs.is_open()) throw runtime_error(file + ".rdm2 cannot be opened");
    string line;
    while (getline(fs, line)) {
      stringstream ss(line);
      int i, j, k, l;
      double re, im;
      // assuming that the 2RDM is dumped as i+ j+ k l -> i k j l
      ss >> i >> j >> k >> l >> re >> im;
      const int ii = ((i-1)/2) + ((i-1)%2)*norb_;
      const int jj = ((j-1)/2) + ((j-1)%2)*norb_;
      const int kk = ((k-1)/2) + ((k-1)%2)*norb_;
      const int ll = ((l-1)/2) + ((l-1)%2)*norb_;
      const complex<double> dat(re, im);
      rdm2_av_expanded_->element(ii, kk, jj, ll) = dat; 
      rdm2_av_expanded_->element(ii, ll, jj, kk) = -dat; 
      rdm2_av_expanded_->element(jj, kk, ii, ll) = -dat; 
      rdm2_av_expanded_->element(jj, ll, ii, kk) = dat; 

      rdm2_av_expanded_->element(kk, ii, ll, jj) = conj(dat);
      rdm2_av_expanded_->element(kk, jj, ll, ii) = -conj(dat);
      rdm2_av_expanded_->element(ll, ii, kk, jj) = -conj(dat);
      rdm2_av_expanded_->element(ll, jj, kk, ii) = conj(dat);
    }
    fs.close();
  }
}
