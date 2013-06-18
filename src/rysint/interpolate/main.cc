//
// Author : Toru Shiozaki
// Date   : May 2009
//

#define NGRID 12
#define MAXT 64
#define NBOX 32
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "mpreal.h"
#include <map>
#include <cmath>
#include <string>
#include <cassert>
#include <fstream>
#include "gmp_macros.h"
#include <boost/lexical_cast.hpp>

extern "C" {
  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
}

using namespace boost;
using namespace std;


void rysroot_gmp(const vector<mpfr::mpreal>& ta, vector<mpfr::mpreal>& dx, vector<mpfr::mpreal>& dw, const int nrank, const int nbatch) ;

using namespace mpfr;
vector<mpreal> chebft(int n) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  vector<mpreal> out(n);
  const mpreal half = "0.5";
  for (int k = 0; k != n; ++k) {
     const mpreal y = mpfr::cos(GMPPI * (k + half) / n);
     out[k] = y;
  }
  return out;
}

vector<vector<double>> get_C(const mpreal tbase, const mpreal stride, int rank) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  const int n = NGRID;

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  vector<mpreal> cheb = chebft(n);

  const mpreal Tmin = tbase;
  const mpreal Tmax = Tmin + stride;
  const mpreal Tp = half * (Tmin + Tmax);

  vector<mpreal> Tpoints(n);
  for (int i = 0; i != n; ++i) {
    Tpoints[i] = stride*half*cheb[i] + Tp;
  }

  vector<map<mpreal, mpreal>> table_reserve(n);
  for (int i = 0; i < n; ++i) {
    vector<mpreal> ttt(1); ttt[0] = Tpoints[i];
    vector<mpreal> dx(rank);
    vector<mpreal> dw(rank);
    rysroot_gmp(ttt, dx, dw, rank, 1);
    // sort dx and dw using dx
    for (int j = 0; j != rank; ++j)
      table_reserve[i].insert(make_pair(dx[j], dw[j]));
  }

  vector<vector<double>> c;
  for (int ii = 0; ii != rank; ++ii) {
    vector<double> tc(n);
    vector<double> tc2(n);

    vector<mpreal> cdx, cdw;
    for (int j = 0; j != n; ++j) {
      auto iter = table_reserve[j].begin();
      for (int i = 0; i != ii; ++i) ++iter;
      cdx.push_back(iter->first);
      cdw.push_back(iter->second);
    }

    const mpreal two = "2.0";
    const mpreal half = "0.5";
    const mpreal fac = two / n;
    const mpreal pi = GMPPI;
    for (int j = 0; j != n; ++j) {
      mpreal sum = "0.0";
      mpreal sum2 = "0.0";
      for (int k = 0; k != n; ++k) {
        sum += cdx[k] * mpfr::cos(pi * j * (k + half) / n);
        sum2 += cdw[k] * mpfr::cos(pi * j * (k + half) / n);
      }
      tc[j] = sum * fac;
      tc2[j] = sum2 * fac;
    }
    if (tc[n-1] > 1.0e-10 || tc2[n-1] > 1.0e-10) {
      cout << " caution: cheb not converged " << ii << " " << setprecision(10) << fixed << (double)Tmin << " " << (double)Tmax << endl;
      for (int i = 0; i != n; ++i) {
        cout << setw(20) << (double)(Tpoints[i]) << setw(20) << (double)cdx[i] << setw(20) << (double)(cdw[i]) << endl;
      }
    }
    c.push_back(tc);
    c.push_back(tc2);
  }
  return c;
}


extern "C" {
  void breitroot1_(double*, double*, double*, const int*);
  void breitroot2_(double*, double*, double*, const int*);
  void breitroot3_(double*, double*, double*, const int*);
  void breitroot4_(double*, double*, double*, const int*);
  void breitroot5_(double*, double*, double*, const int*);
  void breitroot6_(double*, double*, double*, const int*);
  void breitroot7_(double*, double*, double*, const int*);
  void breitroot8_(double*, double*, double*, const int*);
  void breitroot9_(double*, double*, double*, const int*);
  void breitroot10_(double*, double*, double*, const int*);
  void breitroot11_(double*, double*, double*, const int*);
  void breitroot12_(double*, double*, double*, const int*);
  void breitroot13_(double*, double*, double*, const int*);
}

struct B {
  void (*breitroot[14])(double*, double*, double*, const int*);
  B() {
    breitroot[1] = breitroot1_;
    breitroot[2] = breitroot2_;
    breitroot[3] = breitroot3_;
    breitroot[4] = breitroot4_;
    breitroot[5] = breitroot5_;
    breitroot[6] = breitroot6_;
    breitroot[7] = breitroot7_;
    breitroot[8] = breitroot8_;
    breitroot[9] = breitroot9_;
    breitroot[10] = breitroot10_;
    breitroot[11] = breitroot11_;
    breitroot[12] = breitroot12_;
    breitroot[13] = breitroot13_;
  }
  void root(const int rank, double* a, double* b, double* c, const int* d) { breitroot[rank](a,b,c,d); }
} breit;

extern "C" {
  void spin2root1_(double*, double*, double*, const int*);
  void spin2root2_(double*, double*, double*, const int*);
  void spin2root3_(double*, double*, double*, const int*);
  void spin2root4_(double*, double*, double*, const int*);
  void spin2root5_(double*, double*, double*, const int*);
  void spin2root6_(double*, double*, double*, const int*);
  void spin2root7_(double*, double*, double*, const int*);
  void spin2root8_(double*, double*, double*, const int*);
  void spin2root9_(double*, double*, double*, const int*);
  void spin2root10_(double*, double*, double*, const int*);
  void spin2root11_(double*, double*, double*, const int*);
  void spin2root12_(double*, double*, double*, const int*);
  void spin2root13_(double*, double*, double*, const int*);
}

struct S {
  void (*spin2root[14])(double*, double*, double*, const int*);
  S() {
    spin2root[1] = spin2root1_;
    spin2root[2] = spin2root2_;
    spin2root[3] = spin2root3_;
    spin2root[4] = spin2root4_;
    spin2root[5] = spin2root5_;
    spin2root[6] = spin2root6_;
    spin2root[7] = spin2root7_;
    spin2root[8] = spin2root8_;
    spin2root[9] = spin2root9_;
    spin2root[10] = spin2root10_;
    spin2root[11] = spin2root11_;
    spin2root[12] = spin2root12_;
    spin2root[13] = spin2root13_;
  }
  void root(const int rank, double* a, double* b, double* c, const int* d) { spin2root[rank](a,b,c,d); }
} spin2;


bool test(const int nrank, const double tin) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  const static int nsize = 1;
  vector<mpreal> tt(nsize, tin);
  vector<mpreal> rr(nsize*nrank);
  vector<mpreal> ww(nsize*nrank);
  rysroot_gmp(tt, rr, ww, nrank, nsize);
  map<mpreal,mpreal> gmp;
  for (int i = 0; i != nsize*nrank; ++i)
    gmp.insert(make_pair(rr[i], ww[i]));

  double dt[nsize] = {(double)(tt[0])};
  double dr[nsize*nrank];
  double dw[nsize*nrank];
#ifndef SPIN2
#ifdef BREIT
  breit.root(nrank, dt, dr, dw, &nsize);
#else
  assert(false);
#endif
#else
  spin2.root(nrank, dt, dr, dw, &nsize);
#endif

  cout << setprecision(10) << scientific << endl;
  auto iter = gmp.begin();
  for (int i = 0; i != nrank*nsize; ++i, ++iter) {
    cout << setw(20) << dr[i] << setw(20) << iter->first  << setw(20) << fabs(dr[i] - (double)(iter->first)) << endl;
    cout << setw(20) << dw[i] << setw(20) << iter->second << setw(20) << fabs(dw[i] - (double)(iter->second)) << endl;
  }
  iter = gmp.begin();
  for (int i = 0; i != nrank; ++i, ++iter) {
    if (!(fabs(dr[i] - (double)(iter->first)))) cout << dt[0] << endl;
    assert(fabs(dr[i] - (double)(iter->first)) < 1.0e-14);
    assert(fabs(dw[i] - (double)(iter->second)) < 1.0e-14);
  }
  cout << "test passed: rank" << setw(3) << nrank << endl;
}

#include <boost/lexical_cast.hpp>

int main(int argc, char** argv) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  mpfr::mpreal pi = GMPPI;

  if (argc > 1) {
    const string toggle = argv[1];
    if (toggle == "-t") {
#if 0
      if (argc <= 3) assert(false);
      const string low = argv[2];
      const string high = argv[3];
      for (int i = 0; i < 700; ++i) {
        for (int n = lexical_cast<int>(low); n <= lexical_cast<int>(high); ++n) test(n, i*0.1+1.0e-10);
      }
#else
      test(3,1.10033333333333);
      test(3,1.11133333333333);
      test(3,1.12233333333333);
      test(3,1.13233333333333);
      test(3,1.14333333333333);
      test(3,1.14333333333333e3);
#endif
      return 0;
    }
  }

  vector<double> nbox_(14);
  for (int nroot=1; nroot!=14; ++nroot) {
    nbox_[nroot] = NBOX;
  }

  for (int nroot=1; nroot!=14; ++nroot) { // this is the outer most loop.
    if (argc > 2) {
      const string toggle = argv[1];
      if (toggle == "-r") {
        const string target = argv[2];
        if (nroot != lexical_cast<int>(target)) continue;
      }
    }
    vector<double> aroot;
    vector<double> aweight;
#ifndef SPIN2
#ifndef BREIT
    // first obtain asymptotics
    const int n=nroot*2;

    double a[10000] = {0.0};
    double b[100];
    double c[500];
    for (int i=0; i!= n; ++i) {
      a[i+i*n] = 0.0;
      if (i > 0) {
        const double ia = static_cast<double>(i);
        a[(i-1)+i*n] = ::sqrt(ia*0.5);
        a[i+(i-1)*n] = ::sqrt(ia*0.5);
      }
    }
    int nn = n*5;
    int info = 0;
    dsyev_("v", "U", &n, a, &n, b, c, &nn, &info);
    for (int j = 0; j != nroot; ++j) {
      aroot.push_back((double)(b[nroot+j]*b[nroot+j]));
      aweight.push_back((double)(a[(nroot+j)*(nroot*2)]*a[(nroot+j)*(nroot*2)]*sqrt(pi)));
    }
#else
    const mpreal t = 1000;
    const mpreal s = 2000;
    vector<mpreal> dx(nroot*2);
    vector<mpreal> dw(nroot*2);
    vector<mpreal> tt(1, t); tt.push_back(s);
    rysroot_gmp(tt, dx, dw, nroot, 2);
    for (int j = 0; j != nroot; ++j) {
      assert(fabs(dx[j]*t - dx[j+nroot]*s) < 1.0e-16);
      assert(fabs(dw[j]*t*sqrt(t) - dw[j+nroot]*s*sqrt(s)) < 1.0e-16);
      aroot.push_back((double)dx[j]*t);
      aweight.push_back((double)dw[j]*t*sqrt(t));
    }
#endif
#else
    const mpreal t = 1000;
    const mpreal s = 2000;
    vector<mpreal> dx(nroot*2);
    vector<mpreal> dw(nroot*2);
    vector<mpreal> tt(1, t); tt.push_back(s);
    rysroot_gmp(tt, dx, dw, nroot, 2);
    for (int j = 0; j != nroot; ++j) {
      assert(fabs(dx[j]*t - dx[j+nroot]*s) < 1.0e-16);
      assert(fabs(dw[j]*t*t*sqrt(t) - dw[j+nroot]*s*s*sqrt(s)) < 1.0e-16);
      aroot.push_back((double)dx[j]*t);
      aweight.push_back((double)dw[j]*t*t*sqrt(t));
    }
#endif

    const int ndeg = NGRID;
    const int nbox = nbox_[nroot];
    const int jend = nbox;
    const double stride = static_cast<double>(MAXT)/nbox;
    const mpreal mstride = static_cast<mpreal>(MAXT)/nbox;
    ofstream ofs;
#ifndef SPIN2
#ifndef BREIT
    const string func = "eriroot";
#else
    const string func = "breitroot";
#endif
#else
    const string func = "spin2root";
#endif
    string filename = "_" + func + "_" + lexical_cast<string>(nroot) + ".cc";
    ofs.open(filename.c_str());
    ofs << "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) 2013 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
#include <algorithm>" << endl;
#ifndef SPIN2
#ifdef BREIT
ofs << "#include <src/rysint/breitrootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void BreitRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
ofs << "#include <src/rysint/erirootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void ERIRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#endif
#else
ofs << "#include <src/rysint/spin2rootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void Spin2RootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#endif
   ofs << "\
  constexpr double ax["<<nroot<<"] = {";
    for (int j=0; j!= nroot; ++j) {
      ofs << scientific << setprecision(15) << setw(20) << aroot[j];
      if (j != nroot-1) ofs << ",";
      if (j%7 == 4)    ofs << endl << "    ";
    }
    ofs << "};" << endl; 
    ofs << "\
  constexpr double aw["<<nroot<<"] = {";
    for (int j=0; j!= nroot; ++j) {
      ofs << scientific << setprecision(15) << setw(20) << aweight[j];
      if (j != nroot-1) ofs << ",";
      if (j%7 == 4) ofs << endl << "    ";
    }
    ofs << "};" << endl;

////////////////////////////////////////
// now creates data
////////////////////////////////////////
    stringstream listx, listw;
    string indent(" ");
    int nblock = 0;
    int index = 0;
    double tiny = 1.0e-100;
    int xcnt = 0;
    int wcnt = 0;
    for (int j=0; j!= jend; ++j) {
      vector<vector<double>> c_all = get_C(j*mstride,mstride,nroot);

      for (int i = 0; i != nroot; ++i, ++index) {
        const int ii = 2 * i;
        const vector<double> x = c_all[ii];
        const vector<double> w = c_all[ii + 1];

        for (auto iter = x.begin(); iter != x.end(); ++iter) {
          listx << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != x.end() || j+1 != jend || i+1 != nroot) listx << ",";
          if (xcnt++ % 7 == 4) listx << "\n";
        }
        for (auto iter = w.begin(); iter != w.end(); ++iter) {
          listw << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != w.end() || j+1 != jend || i+1 != nroot) listw << ",";
          if (wcnt++ % 7 == 4) listw << "\n";
        }
      }
    }
#ifndef SPIN2
#ifndef BREIT
    string tafactor = "t";
#else
    string tafactor = "t*t*t";
#endif
#else
    string tafactor = "t*t*t*t*t";
#endif
    ofs << "\
  constexpr double x[" << nroot*nbox*ndeg<<"] = {";
    ofs << listx.str() << "\
  };" << endl;
    ofs << "\
  constexpr double w[" << nroot*nbox*ndeg<<"] = {";
    ofs << listw.str() << "\
  };" << endl;

  ofs << "\
  int offset = -" << nroot << ";\n\
  for (int i = 1; i <= n; ++i) {\n\
    double t = ta[i-1];\n\
    offset += " << nroot << ";\n\
    if (t < 0.0) {\n\
      fill_n(rr+offset, " << nroot << ", 0.5);\n\
      fill_n(ww+offset, " << nroot << ", 0.0);\n\
    } else if (t >= " << MAXT << ".0) {\n\
      t = 1.0/sqrt(t);\n\
      for (int r = 0; r != " << nroot << "; ++r) {\n\
        rr[offset+r] = ax[r]*t*t;\n\
        ww[offset+r] = aw[r]*" + tafactor + ";\n\
      }\n\
    } else {\n\
      int it = static_cast<int>(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< ");\n\
      t = (t-it*" << stride << "-" << setw(20) << setprecision(15) << fixed << stride/2.0 << ") *" << setw(20) << setprecision(15) << fixed << 2.0/stride << ";\n\
      const double t2 = t * 2.0;\n\
      for (int j=1; j <=" << nroot << "; ++j) {\n\
        const int boxof = it*" << ndeg*nroot << "+" << ndeg << "*(j-1);\n";
     assert((ndeg/2)*2 == ndeg);
     for (int i=ndeg; i!=0; --i) {
       if (i==ndeg) {
         ofs << "\
        double d = x[boxof+" << i-1 << "];\n\
        double e = w[boxof+" << i-1 << "];\n";
       } else if (i==ndeg-1) {
         ofs << "\
        double f = t2*d + x[boxof+" << i-1 << "];\n\
        double g = t2*e + w[boxof+" << i-1 << "];\n";
       } else if (i != 1 && ((i/2)*2 == i)) { // odd
         ofs << "\
        d = t2*f - d + x[boxof+" << i-1 << "];\n\
        e = t2*g - e + w[boxof+" << i-1 << "];\n";
       } else if (i != 1) { // even
         ofs << "\
        f = t2*d - f + x[boxof+" << i-1 << "];\n\
        g = t2*e - g + w[boxof+" << i-1 << "];\n";
       } else {
         ofs << "\
        rr[offset+j-1] = t*d - f + x[boxof+" << i-1 << "]*0.5;\n\
        ww[offset+j-1] = t*e - g + w[boxof+" << i-1 << "]*0.5;\n";
       }
     }

     ofs << "\
      }\n\
    }\n\
  }\n\
}";

    ofs.close();
  }


  return 0;
}

