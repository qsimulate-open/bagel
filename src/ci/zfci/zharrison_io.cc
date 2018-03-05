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
      const int jfac  = (i & 1) + 1;
      const int j2fac = (i & 2)/2 + 1;
      const int kfac  = (i & 4)/4 + 1;
      const int k2fac = (i & 8)/8 + 1;

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

    for (int i = 0; i != 4; ++i) {
      if (!jop_->mo1e()->exist(i)) continue;
      cout << "Writing 1e integral block " << i+1 << " / 4" << endl;
      shared_ptr<const ZMatrix> tmp = jop_->mo1e(i);

      const int jfac = (i & 1) + 1;
      const int kfac = (i & 2)/2 + 1;

      for (int j = 0; j != norb_; ++j)
        for (int k = 0; k != norb_; ++k) {
          const complex<double> val = tmp->element(k, j);
          if (abs(val) > 1.0e-9)
            fs << val << setw(4) << 2*k+kfac << setw(4) << 2*j+jfac << "   0   0" << endl;
        }
    }
    fs << "(" << jop_->core_energy()  + geom_->nuclear_repulsion() << ",0.0)" << "   0   0   0   0" << endl;
    fs.close();
  }
  else throw runtime_error("Unable to open file: ZHarrison::dump_ints");
}


std::shared_ptr<Kramers<8,ZRDM<4>>> ZHarrison::read_external_rdm4(const int ist, const int jst, const string& file) const {
  auto out = make_shared<Kramers<8,ZRDM<4>>>();

  map<array<int,4>,double> elem;
  elem.emplace(array<int,4>{{0,1,2,3}},  1.0); elem.emplace(array<int,4>{{0,1,3,2}}, -1.0); elem.emplace(array<int,4>{{0,2,1,3}}, -1.0);
  elem.emplace(array<int,4>{{0,2,3,1}},  1.0); elem.emplace(array<int,4>{{0,3,1,2}},  1.0); elem.emplace(array<int,4>{{0,3,2,1}}, -1.0);
  elem.emplace(array<int,4>{{1,0,2,3}}, -1.0); elem.emplace(array<int,4>{{1,0,3,2}},  1.0); elem.emplace(array<int,4>{{1,2,0,3}},  1.0);
  elem.emplace(array<int,4>{{1,2,3,0}}, -1.0); elem.emplace(array<int,4>{{1,3,0,2}}, -1.0); elem.emplace(array<int,4>{{1,3,2,0}},  1.0);
  elem.emplace(array<int,4>{{2,0,1,3}},  1.0); elem.emplace(array<int,4>{{2,0,3,1}}, -1.0); elem.emplace(array<int,4>{{2,1,0,3}}, -1.0);
  elem.emplace(array<int,4>{{2,1,3,0}},  1.0); elem.emplace(array<int,4>{{2,3,0,1}},  1.0); elem.emplace(array<int,4>{{2,3,1,0}}, -1.0);
  elem.emplace(array<int,4>{{3,0,1,2}}, -1.0); elem.emplace(array<int,4>{{3,0,2,1}},  1.0); elem.emplace(array<int,4>{{3,1,0,2}},  1.0);
  elem.emplace(array<int,4>{{3,1,2,0}}, -1.0); elem.emplace(array<int,4>{{3,2,0,1}}, -1.0); elem.emplace(array<int,4>{{3,2,1,0}},  1.0);

  stringstream ss; ss << file << "_" << ist << "_" << jst << ".rdm4";
  ifstream fs(ss.str());
  if (!fs.is_open()) throw runtime_error(ss.str() + " cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j, k, l, m, n, o, p;
    double re, im;
    // assuming that the 2RDM is dumped as i+ j+ k+ l m n -> i l j m k n
    ss >> i >> j >> k >> o >> l >> m >> n >> p >> re >> im;
    assert(i <= norb_*2 && j <= norb_*2 && k <= norb_*2 && l <= norb_*2 && m <= norb_*2 && n <= norb_*2 && o <= norb_*2 && p <= norb_*2);
    map<int,pair<int,int>> mij{{0,{(i-1)/2,(i-1)%2}}, {1,{(j-1)/2,(j-1)%2}}, {2,{(k-1)/2,(k-1)%2}}, {3,{(o-1)/2,(o-1)%2}}};
    map<int,pair<int,int>> mkl{{0,{(l-1)/2,(l-1)%2}}, {1,{(m-1)/2,(m-1)%2}}, {2,{(n-1)/2,(n-1)%2}}, {3,{(p-1)/2,(p-1)%2}}};
    const complex<double> dat(re, im);
    for (auto& eij : elem) {
      for (auto& ekl : elem) {
        if (mij[eij.first[0]].second > mij[eij.first[1]].second || mij[eij.first[1]].second > mij[eij.first[2]].second || mij[eij.first[2]].second > mij[eij.first[3]].second ||
            mkl[ekl.first[0]].second > mkl[ekl.first[1]].second || mkl[ekl.first[1]].second > mkl[ekl.first[2]].second || mkl[ekl.first[2]].second > mkl[ekl.first[3]].second) continue;

        const KTag<8> t{mij[eij.first[0]].second, mkl[ekl.first[0]].second, mij[eij.first[1]].second, mkl[ekl.first[1]].second,
                        mij[eij.first[2]].second, mkl[ekl.first[2]].second, mij[eij.first[3]].second, mkl[ekl.first[3]].second};
        if (!out->exist(t))
          out->add(t, make_shared<ZRDM<4>>(norb_));
        out->at(t)->element(mij[eij.first[0]].first, mkl[ekl.first[0]].first, mij[eij.first[1]].first, mkl[ekl.first[1]].first,
                            mij[eij.first[2]].first, mkl[ekl.first[2]].first, mij[eij.first[3]].first, mkl[ekl.first[3]].first)
          = eij.second * ekl.second * dat;

        if (ist == jst) {
          const KTag<8> t2{mkl[ekl.first[0]].second, mij[eij.first[0]].second, mkl[ekl.first[1]].second, mij[eij.first[1]].second,
                           mkl[ekl.first[2]].second, mij[eij.first[2]].second, mkl[ekl.first[3]].second, mij[eij.first[3]].second};
          if (!out->exist(t2))
            out->add(t2, make_shared<ZRDM<4>>(norb_));
          out->at(t2)->element(mkl[ekl.first[0]].first, mij[eij.first[0]].first, mkl[ekl.first[1]].first, mij[eij.first[1]].first,
                               mkl[ekl.first[2]].first, mij[eij.first[2]].first, mkl[ekl.first[3]].first, mij[eij.first[3]].first)
            = eij.second * ekl.second * conj(dat);
        }
      }
    }
  }
  for (auto& i : elem) {
    for (auto& j : elem) {
      vector<int> perm(8);
      for (int k = 0; k != 4; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      out->emplace_perm(perm, j.second*i.second);
    }
  }
  return out;
}


std::shared_ptr<Kramers<6,ZRDM<3>>> ZHarrison::read_external_rdm3(const int ist, const int jst, const string& file, const bool fock_contracted) const {
  auto out = make_shared<Kramers<6,ZRDM<3>>>();

  map<array<int,3>,double> elem;
  elem.emplace(array<int,3>{{0,1,2}},  1.0); elem.emplace(array<int,3>{{0,2,1}}, -1.0); elem.emplace(array<int,3>{{1,0,2}}, -1.0);
  elem.emplace(array<int,3>{{1,2,0}},  1.0); elem.emplace(array<int,3>{{2,0,1}},  1.0); elem.emplace(array<int,3>{{2,1,0}}, -1.0);

  stringstream ss; ss << file << "_" << ist << "_" << jst << (fock_contracted ? ".rdm4f" : ".rdm3");
  ifstream fs(ss.str());
  if (!fs.is_open()) throw runtime_error(ss.str() + " cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j, k, l, m, n;
    double re, im;
    // assuming that the 2RDM is dumped as i+ j+ k+ l m n -> i l j m k n
    ss >> i >> j >> k >> l >> m >> n >> re >> im;
    assert(i <= norb_*2 && j <= norb_*2 && k <= norb_*2 && l <= norb_*2 && m <= norb_*2 && n <= norb_*2);
    map<int,pair<int,int>> mij{{0,{(i-1)/2,(i-1)%2}}, {1,{(j-1)/2,(j-1)%2}}, {2,{(k-1)/2,(k-1)%2}}};
    map<int,pair<int,int>> mkl{{0,{(l-1)/2,(l-1)%2}}, {1,{(m-1)/2,(m-1)%2}}, {2,{(n-1)/2,(n-1)%2}}};
    const complex<double> dat(re, im);
    for (auto& eij : elem) {
      for (auto& ekl : elem) {
        if (mij[eij.first[0]].second > mij[eij.first[1]].second || mij[eij.first[1]].second > mij[eij.first[2]].second ||
            mkl[ekl.first[0]].second > mkl[ekl.first[1]].second || mkl[ekl.first[1]].second > mkl[ekl.first[2]].second) continue;

        const KTag<6> t{mij[eij.first[0]].second, mkl[ekl.first[0]].second, mij[eij.first[1]].second, mkl[ekl.first[1]].second,
                        mij[eij.first[2]].second, mkl[ekl.first[2]].second};
        if (!out->exist(t))
          out->add(t, make_shared<ZRDM<3>>(norb_));
        out->at(t)->element(mij[eij.first[0]].first, mkl[ekl.first[0]].first, mij[eij.first[1]].first, mkl[ekl.first[1]].first,
                            mij[eij.first[2]].first, mkl[ekl.first[2]].first)
          = eij.second * ekl.second * dat;

        if (ist == jst) {
          const KTag<6> t2{mkl[ekl.first[0]].second, mij[eij.first[0]].second, mkl[ekl.first[1]].second, mij[eij.first[1]].second,
                           mkl[ekl.first[2]].second, mij[eij.first[2]].second};
          if (!out->exist(t2))
            out->add(t2, make_shared<ZRDM<3>>(norb_));
          out->at(t2)->element(mkl[ekl.first[0]].first, mij[eij.first[0]].first, mkl[ekl.first[1]].first, mij[eij.first[1]].first,
                               mkl[ekl.first[2]].first, mij[eij.first[2]].first)
            = eij.second * ekl.second * conj(dat);
        }
      }
    }
  }
  for (auto& i : elem) {
    for (auto& j : elem) {
      vector<int> perm(6);
      for (int k = 0; k != 3; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      out->emplace_perm(perm, j.second*i.second);
    }
  }
  return out;
}


std::shared_ptr<Kramers<4,ZRDM<2>>> ZHarrison::read_external_rdm2(const int ist, const int jst, const string& file) const {
  auto out = make_shared<Kramers<4,ZRDM<2>>>();

  map<array<int,2>,double> elem;
  elem.emplace(array<int,2>{{0,1}},  1.0); elem.emplace(array<int,2>{{1,0}}, -1.0);

  stringstream ss; ss << file << "_" << ist << "_" << jst << ".rdm2";
  ifstream fs(ss.str());
  if (!fs.is_open()) throw runtime_error(ss.str() + " cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j, k, l;
    double re, im;
    // assuming that the 2RDM is dumped as i+ j+ k l -> i k j l
    ss >> i >> j >> k >> l >> re >> im;
    assert(i <= norb_*2 && j <= norb_*2 && k <= norb_*2 && l <= norb_*2);
    map<int,pair<int,int>> mij{{0,{(i-1)/2,(i-1)%2}}, {1,{(j-1)/2,(j-1)%2}}};
    map<int,pair<int,int>> mkl{{0,{(k-1)/2,(k-1)%2}}, {1,{(l-1)/2,(l-1)%2}}};
    const complex<double> dat(re, im);
    for (auto& eij : elem) {
      for (auto& ekl : elem) {
        if (mij[eij.first[0]].second > mij[eij.first[1]].second || mkl[ekl.first[0]].second > mkl[ekl.first[1]].second) continue;

        const KTag<4> t{mij[eij.first[0]].second, mkl[ekl.first[0]].second, mij[eij.first[1]].second, mkl[ekl.first[1]].second};
        if (!out->exist(t))
          out->add(t, make_shared<ZRDM<2>>(norb_));
        out->at(t)->element(mij[eij.first[0]].first, mkl[ekl.first[0]].first, mij[eij.first[1]].first, mkl[ekl.first[1]].first)
          = eij.second * ekl.second * dat;

        if (ist == jst) {
          const KTag<4> t2{mkl[ekl.first[0]].second, mij[eij.first[0]].second, mkl[ekl.first[1]].second, mij[eij.first[1]].second};
          if (!out->exist(t2))
            out->add(t2, make_shared<ZRDM<2>>(norb_));
          out->at(t2)->element(mkl[ekl.first[0]].first, mij[eij.first[0]].first, mkl[ekl.first[1]].first, mij[eij.first[1]].first)
            = eij.second * ekl.second * conj(dat);
        }
      }
    }
  }
  for (auto& i : elem) {
    for (auto& j : elem) {
      vector<int> perm(4);
      for (int k = 0; k != 2; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      out->emplace_perm(perm, j.second*i.second);
    }
  }
  return out;
}

std::shared_ptr<Kramers<2,ZRDM<1>>> ZHarrison::read_external_rdm1(const int ist, const int jst, const string& file) const {
  auto out = make_shared<Kramers<2,ZRDM<1>>>();
  auto tmp = make_shared<ZRDM<1>>(norb_);
  out->add(KTag<2>("00"), tmp->clone());
  out->add(KTag<2>("01"), tmp->clone());
  out->add(KTag<2>("10"), tmp->clone());
  out->add(KTag<2>("11"), tmp);

  stringstream ss; ss << file << "_" << ist << "_" << jst << ".rdm1";
  ifstream fs(ss.str());
  if (!fs.is_open()) throw runtime_error(ss.str() + " cannot be opened");
  string line;
  while (getline(fs, line)) {
    stringstream ss(line);
    int i, j;
    double re, im;
    ss >> i >> j >> re >> im;
    assert(i <= norb_*2 && j <= norb_*2);

    const int ii = (i-1)/2;
    const int jj = (j-1)/2;
    const complex<double> dat(re, im);
    const int ti = (i-1)%2;
    const int tj = (j-1)%2;
    const KTag<2> tag1{ti, tj};
    out->at(tag1)->element(ii, jj) = dat;
    if (ist == jst) {
      const KTag<2> tag2{tj, ti};
      out->at(tag2)->element(jj, ii) = conj(dat);
    }
  }
  return out;
}


void ZHarrison::read_external_rdm12_av(const string& file) {
  // feed RDM1
  rdm1_av_ = read_external_rdm1(/*TODO*/0, 0, file);
  rdm1_av_expanded_ = expand_kramers<1,complex<double>>(rdm1_av_, norb_);
  rdm2_av_ = read_external_rdm2(/*TODO*/0, 0, file);
  rdm2_av_expanded_ = expand_kramers<2,complex<double>>(rdm2_av_, norb_);
}
