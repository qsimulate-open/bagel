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

#include "../r2rootlist.h"

extern "C" {
  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
}

using namespace boost;
using namespace std;
using namespace mpfr;
using namespace bagel;

void rysroot_gmp(const vector<mpfr::mpreal>& ta, vector<mpfr::mpreal>& dx, vector<mpfr::mpreal>& dw, const int nrank, const int nbatch) ;

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

vector<vector<double>> get_C(const mpreal tbase, const mpreal stride, int rank, const bool asymp) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  const int n = NGRID;

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const mpreal one = "1.0";
  vector<mpreal> cheb = chebft(n);

  const mpreal Tmin = tbase;
  const mpreal Tmax = Tmin + stride;
  const mpreal Tp = half * (Tmin + Tmax);

  vector<mpreal> Tpoints(n);
  for (int i = 0; i != n; ++i) {
    Tpoints[i] = stride*half*cheb[i] + Tp;
  }
  if (asymp) {
    for (int i = 0; i != n; ++i) {
      Tpoints[i] = one / (Tpoints[i] * Tpoints[i]);
    }
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
      tc[j] = (sum * fac).toDouble();
      tc2[j] = (sum2 * fac).toDouble();
    }
    if (tc[n-1] > 1.0e-10 || tc2[n-1] > 1.0e-10) {
      cout << " caution: cheb not converged " << ii << " " << setprecision(10) << fixed << Tmin.toDouble() << " " << Tmax.toDouble() << endl;
      for (int i = 0; i != n; ++i) {
        //cout << setw(20) << Tpoints[i].toDouble() << setw(20) << cdx[i].toDouble() << setw(20) << cdw[i].toDouble() << endl;
        cout << setw(20) << Tpoints[i].toDouble() << setw(20) << tc[n-1] << setw(20) << tc2[n-1] << endl;
      }
    }
    c.push_back(tc);
    c.push_back(tc2);
  }
  return c;
}


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

  double dt[nsize] = {tt[0].toDouble()};
  double dr[nsize*nrank];
  double dw[nsize*nrank];

  r2root__.root(nrank, dt, dr, dw, nsize);

  cout << setprecision(10) << scientific << endl;
  auto iter = gmp.begin();
  for (int i = 0; i != nrank*nsize; ++i, ++iter) {
    cout << setw(20) << dr[i] << setw(20) << iter->first  << setw(20) << fabs(dr[i] - (iter->first).toDouble()) << endl;
    cout << setw(20) << dw[i] << setw(20) << iter->second << setw(20) << fabs(dw[i] - (iter->second).toDouble()) << endl;
  }
  iter = gmp.begin();
  for (int i = 0; i != nrank; ++i, ++iter) {
    if (!(fabs(dr[i] - (iter->first).toDouble()))) cout << dt[0] << endl;
    assert(fabs(dr[i] - (iter->first).toDouble()) < 1.0e-13);
    assert(fabs(dw[i] - (iter->second).toDouble()) < 1.0e-13);
  }
  cout << "test passed: rank" << setw(3) << nrank << endl;
}

#include <boost/lexical_cast.hpp>

int main(int argc, char** argv) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  mpfr::mpreal pi = GMPPI;

  if (argc > 1) {
    cout << "--- TEST---" << endl;
    const string toggle = argv[1];
    if (toggle == "t") {
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
      test(3,1.14333333333333e1);
      test(3,1.14333333333333e3);
      test(3,1.14333333333333e4);
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
#ifndef DAWSON
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
        a[(i-1)+i*n] = std::sqrt(ia*0.5);
        a[i+(i-1)*n] = std::sqrt(ia*0.5);
      }
    }
    int nn = n*5;
    int info = 0;
    dsyev_("v", "U", &n, a, &n, b, c, &nn, &info);
    for (int j = 0; j != nroot; ++j) {
      aroot.push_back(b[nroot+j]*b[nroot+j]);
      aweight.push_back(a[(nroot+j)*(nroot*2)]*a[(nroot+j)*(nroot*2)]*(sqrt(pi)).toDouble());
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
      aroot.push_back((dx[j]*t).toDouble());
      aweight.push_back((dw[j]*t*sqrt(t)).toDouble());
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
      aroot.push_back((dx[j]*t).toDouble());
      aweight.push_back((dw[j]*t*t*sqrt(t)).toDouble());
    }
#endif
#endif

    const int ndeg = NGRID;
    const int nbox = nbox_[nroot];
#ifndef DAWSON
    const int jend = nbox;
#else
    const int jend = nbox + 1;
#endif
    const double stride = static_cast<double>(MAXT)/nbox;
    const mpreal mstride = static_cast<mpreal>(MAXT)/nbox;
    ofstream ofs;
#ifndef SPIN2
#ifndef BREIT
#ifndef DAWSON
    const string func = "eriroot";
#else
    const string func = "r2root";
#endif
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
// the Free Software Foundation; either version 3, or (at your option)\n\
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
ofs << "#include <src/integral/rys/breitrootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void BreitRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
#ifndef DAWSON
ofs << "#include <src/integral/rys/erirootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void ERIRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
ofs << "#include <src/integral/rys/r2rootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void R2RootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#endif
#endif
#else
ofs << "#include <src/integral/rys/spin2rootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void Spin2RootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#endif
#ifndef DAWSON
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
#endif

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
    for (int j=0; j != jend; ++j) {
#ifndef DAWSON
      vector<vector<double>> c_all = get_C(j*mstride, mstride, nroot, false);
#else
      vector<vector<double>> c_all;
      if (j != jend-1) {
        c_all = get_C(j*mstride, mstride, nroot, false);
      } else {
        const mpreal zero = "0.0";
        const mpreal one = "1.0";
        c_all = get_C(zero, one / sqrt(MAXT), nroot, true);
      }
#endif

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

#ifndef DAWSON
    const int nbox1 = nbox;
#else
    const int nbox1 = nbox + 1;
#endif
  ofs << "\
  constexpr double x[" << nroot*nbox1*ndeg<<"] = {";
  ofs << listx.str() << "\
  };" << endl;
  ofs << "\
  constexpr double w[" << nroot*nbox1*ndeg<<"] = {";
  ofs << listw.str() << "\
  };" << endl;

  ofs << "\
  int offset = -" << nroot << ";\n\
  for (int i = 1; i <= n; ++i) {\n\
    double t = ta[i-1];\n\
    offset += " << nroot << ";\n\
    if (t < 0.0) {\n\
      fill_n(rr+offset, " << nroot << ", 0.5);\n\
      fill_n(ww+offset, " << nroot << ", 0.0);\n";
#ifndef DAWSON
      ofs << "\
    } else if (t >= " << MAXT << ".0) {\n\
      t = 1.0/sqrt(t);\n\
      for (int r = 0; r != " << nroot << "; ++r) {\n\
        rr[offset+r] = ax[r]*t*t;\n\
        ww[offset+r] = aw[r]*" + tafactor + ";\n\
      }\n";
#endif
      ofs << "\
    } else {\n";
#ifndef DAWSON
      ofs << "\
      int it = static_cast<int>(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< ");\n";
#else
      ofs << "\
      int it; \n\
      if (t >= " << MAXT << ".0) { \n\
        t = " << setw(20) << setprecision(15) << fixed << 1.0 << " / (t * t) ; \n\
        it = static_cast<int>(" << MAXT << "*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< " + " << setw(20) << setprecision(15) << fixed << stride/2.0 << ");\n\
      } else { \n\
        it = static_cast<int>(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< ");\n\
      } \n";
#endif
      ofs << "\
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

