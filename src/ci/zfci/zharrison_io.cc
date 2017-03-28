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

#include <map>
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
      // assuming here that the fastest bit in i corresponds to the slowest orbital in mo2e
      const int jfac  = 2 - (i & 1);
      const int j2fac = 2 - (i & 2)/2;
      const int kfac  = 2 - (i & 4)/4;
      const int k2fac = 2 - (i & 8)/8;

      const map<int,string> pm{{1, "+"}, {2, "-"}};
      cout << "( " << pm.at(k2fac) << " " << pm.at(j2fac) << " | " << pm.at(kfac) << " " << pm.at(jfac) << " )" << endl;

      for (int j = 0; j != norb_; ++j)
        for (int j2 = 0; j2 != norb_; ++j2)
          for (int k = 0; k != norb_; ++k)
            for (int k2 = 0; k2 != norb_; ++k2) {
              const complex<double> val = tmp->element(k2+norb_*k, j2+norb_*j);
              if (abs(val) > 1.0e-9) {
                fs << setw(20) << val << setw(4) << 2*k2+k2fac << setw(4)
                  << 2*j2+j2fac << setw(4) << 2*k+kfac << setw(4) << 2*j+jfac << endl;   // electron 1, electron 2
              }
            }
    }

    complex<double> tval = 0.0;
    for (int i = 0; i != 4; ++i) {
      if (!jop_->mo1e()->exist(i)) continue;
      cout << "Writing 1e integral block " << i+1 << " / 4" << endl;
      shared_ptr<const ZMatrix> tmp = jop_->mo1e(i);

      const int jfac = 2 - (i & 1);
      const int kfac = 2 - (i & 2)/2;

      for (int j = 0; j != norb_; ++j)
        for (int k = 0; k != norb_; ++k) {
          const complex<double> val = tmp->element(k, j);
          if (abs(val) > 1.0e-9)
            fs << val << setw(4) << 2*k+kfac << setw(4) << 2*j+jfac << "   0   0" << endl;
          if ((j == k) && ((i == 0) || (i == 3))) tval += val;
        }
    }
    fs << "(" << jop_->core_energy()  + geom_->nuclear_repulsion() << ",0.0)" << "   0   0   0   0" << endl;
    fs.close();
  }
  else throw runtime_error("Unable to open file: ZHarrison::dump_ints");
}


std::shared_ptr<Kramers<4,ZRDM<2>>> ZHarrison::read_external_rdm2(const string& file) const {
  auto out = make_shared<Kramers<4,ZRDM<2>>>();
  auto tmp = make_shared<ZRDM<2>>(norb_); 
  out->add(KTag<4>("0000"), tmp->clone());
  out->add(KTag<4>("0001"), tmp->clone());
  out->add(KTag<4>("0101"), tmp->clone());
  out->add(KTag<4>("0010"), tmp->clone());
  out->add(KTag<4>("0011"), tmp->clone());
  out->add(KTag<4>("0111"), tmp->clone());
  out->add(KTag<4>("1010"), tmp->clone());
  out->add(KTag<4>("1011"), tmp->clone());
  out->add(KTag<4>("1111"), tmp);

  ifstream fs(file + ".rdm2"); 
  if (!fs.is_open()) throw runtime_error(file + ".rdm2 cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j, k, l;
    double re, im;
    // assuming that the 2RDM is dumped as i+ j+ k l -> i k j l
    ss >> i >> j >> k >> l >> re >> im;
    const int ii = (i-1)/2;
    const int jj = (j-1)/2;
    const int kk = (k-1)/2;
    const int ll = (l-1)/2;
    const int ti = (i-1)%2;
    const int tj = (j-1)%2;
    const int tk = (k-1)%2;
    const int tl = (l-1)%2;
    const complex<double> dat(re, im);
    {
      const KTag<4> t{ti, tk, tj, tl}; 
      if (out->exist(t)) out->at(t)->element(ii, kk, jj, ll) = dat; 
    } {
      const KTag<4> t{ti, tl, tj, tk}; 
      if (out->exist(t)) out->at(t)->element(ii, ll, jj, kk) = -dat; 
    } {
      const KTag<4> t{tj, tk, ti, tl};
      if (out->exist(t)) out->at(t)->element(jj, kk, ii, ll) = -dat; 
    } {
      const KTag<4> t{tj, tl, ti, tk};
      if (out->exist(t)) out->at(t)->element(jj, ll, ii, kk) = dat; 
    } {
      const KTag<4> t{tk, ti, tl, tj};
      if (out->exist(t)) out->at(t)->element(kk, ii, ll, jj) = conj(dat);
    } {
      const KTag<4> t{tk, tj, tl, ti};
      if (out->exist(t)) out->at(t)->element(kk, jj, ll, ii) = -conj(dat);
    } {
      const KTag<4> t{tl, ti, tk, tj};
      if (out->exist(t)) out->at(t)->element(ll, ii, kk, jj) = -conj(dat);
    } {
      const KTag<4> t{tl, tj, tk, ti};
      if (out->exist(t)) out->at(t)->element(ll, jj, kk, ii) = conj(dat);
    }
  }

  out->emplace_perm({0,1,2,3},  1.0);
  out->emplace_perm({0,3,2,1}, -1.0);
  out->emplace_perm({2,1,0,3}, -1.0);
  out->emplace_perm({2,3,0,1},  1.0);
  return out;
}

std::shared_ptr<Kramers<2,ZRDM<1>>> ZHarrison::read_external_rdm1(const string& file) const {
  auto out = make_shared<Kramers<2,ZRDM<1>>>();
  auto tmp = make_shared<ZRDM<1>>(norb_); 
  out->add(KTag<2>("00"), tmp->clone());
  out->add(KTag<2>("01"), tmp->clone());
  out->add(KTag<2>("10"), tmp->clone());
  out->add(KTag<2>("11"), tmp);

  ifstream fs(file + ".rdm1"); 
  if (!fs.is_open()) throw runtime_error(file + ".rdm1 cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j;
    double re, im;
    ss >> i >> j >> re >> im;

    const int ii = (i-1)/2;
    const int jj = (j-1)/2;
    const complex<double> dat(re, im);
    const int ti = (i-1)%2;
    const int tj = (j-1)%2;
    const KTag<2> tag1{ti, tj};
    const KTag<2> tag2{tj, ti};
    assert(out->exist(tag1));
    assert(out->exist(tag2));
    out->at(tag1)->element(ii, jj) = dat; 
    out->at(tag2)->element(jj, ii) = conj(dat);
  }
  return out;
}


void ZHarrison::read_external_rdm12_av(const string& file) {
  rdm2_av_expanded_ = make_shared<ZRDM<2>>(norb_*2);

  // feed RDM1
  rdm1_av_ = read_external_rdm1(file);
  rdm1_av_expanded_ = expand_kramers<1,complex<double>>(rdm1_av_, norb_);
  rdm2_av_ = read_external_rdm2(file);
  rdm2_av_expanded_ = expand_kramers<2,complex<double>>(rdm2_av_, norb_);
}
