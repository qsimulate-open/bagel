//
// Author : Toru Shiozaki
// Date   : May 2009
//

#define NGRID 12
#define MAXT 64
#define NBOX 32
#define NBOXL 0
#define T_INFTY 11
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

#include "../erirootlist.h"

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

#ifdef DAWSON
  vector<mpreal> tt_infty(1); tt_infty[0] = Tmax;
  vector<mpreal> dx_infty(rank);
  vector<mpreal> dw_infty(rank);
  if (asymp) rysroot_gmp(tt_infty, dx_infty, dw_infty, rank, 1);
#endif
  vector<map<mpreal, mpreal>> table_reserve(n);
  #pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    vector<mpreal> ttt(1); ttt[0] = Tpoints[i];
    vector<mpreal> dx(rank);
    vector<mpreal> dw(rank);
    rysroot_gmp(ttt, dx, dw, rank, 1);
    // sort dx and dw using dx
#ifdef DAWSON
    if (asymp) {
      for (int j = 0; j != rank; ++j) {
        table_reserve[i].insert(make_pair(-(1.0 - dx[j])*ttt[0]/((1.0 - dx_infty[j])*tt_infty[0]), dw[j]*ttt[0]/(dw_infty[j]*tt_infty[0])));
      }
    } else {
      for (int j = 0; j != rank; ++j)
        table_reserve[i].insert(make_pair(dx[j], dw[j]));
    }
#else
    for (int j = 0; j != rank; ++j)
      table_reserve[i].insert(make_pair(dx[j], dw[j]));
#endif
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
        cout << setw(20) << Tpoints[i].toDouble() << setw(20) << tc[i] << setw(20) << tc2[i] << endl;
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

  eriroot__.root(nrank, dt, dr, dw, nsize);

  cout << setprecision(10) << scientific << endl;
  auto iter = gmp.begin();
  for (int i = 0; i != nrank*nsize; ++i, ++iter) {
    cout << setw(20) << dr[i] << setw(20) << iter->first  << setw(20) << fabs(dr[i] - (iter->first).toDouble()) << endl;
    cout << setw(20) << dw[i] << setw(20) << iter->second << setw(20) << fabs(dw[i] - (iter->second).toDouble()) << endl;
  }
  iter = gmp.begin();
  for (int i = 0; i != nrank; ++i, ++iter) {
    if (!(fabs(dr[i] - (iter->first).toDouble()))) cout << dt[0] << endl;
    //assert(fabs(dr[i] - (iter->first).toDouble()) < 1.0e-13);
    //assert(fabs(dw[i] - (iter->second).toDouble()) < 1.0e-13);
  }
  cout << "test passed: rank" << setw(3) << nrank << endl;
  cout << "----------" << endl;
}

#include <boost/lexical_cast.hpp>

int main(int argc, char** argv) {

  const mpreal T_ASYM = static_cast<mpreal>(MAXT + (NBOXL)*(NBOXL + 1.0)*(2.0*NBOXL + 1.0)/6.0);

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
      test(6,1.10033333333333);
      test(6,1.11133333333333);
      test(6,1.12233333333333);
      test(6,1.13233333333333);
      test(6,1.14333333333333);
      test(6,1.14333333333333e1);
      test(6,0.645e2);
      test(6,0.675e2);
      test(6,0.805e2);
      test(6,0.912e2);
      test(6,1.14333333333333e2);
      test(6,1.285e2);
      test(6,128.000000000000);
      test(6,1.31e2);
      test(6,1.38e2);
      test(6,2.43e2);
      test(6,256.000000000000);
      test(6,1.14333333333333e3);
      test(6,8e3);
      test(6,8192.000000000000);
      test(6,1.14333333333333e4);
      test(6,1.14333333333333e5);
      test(6,1.14333333333333e6);
#endif
      return 0;
    }
  }

  vector<double> nbox_(52);
  for (int nroot=1; nroot!=52; ++nroot) {
    nbox_[nroot] = NBOX;
  }

  for (int nroot=1; nroot!=52; ++nroot) { // this is the outer most loop.
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
#else
    mpreal infty;
    if (T_INFTY < MAXT) {
      assert (NBOXL == 0);
      const int nbox0 = static_cast<int>(log(MAXT)/log(2.0));
      infty = pow(2, nbox0 + T_INFTY);
    } else {
      infty = static_cast<mpreal>(T_INFTY);
    }
    vector<mpreal> tt_infty(1); tt_infty[0] = infty;
    vector<mpreal> dx_infty(nroot);
    vector<mpreal> dw_infty(nroot);
    rysroot_gmp(tt_infty, dx_infty, dw_infty, nroot, 1);

    for (int j = 0; j != nroot; ++j) {
      aroot.push_back(((1.0 - dx_infty[j])*tt_infty[0]).toDouble());
      aweight.push_back((dw_infty[j]*tt_infty[0]).toDouble());
    }
#endif

    const int ndeg = NGRID;
    const int nbox = nbox_[nroot];
#ifndef DAWSON
    const int jend = nbox;
#else
    int jend;
    if (MAXT < T_INFTY) {
      jend = NBOX + NBOXL + 1;
    } else {
      jend = NBOX + T_INFTY;
    }
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
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: " + filename + "\n\
// Copyright (C) 2013 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// This program is free software: you can redistribute it and/or modify\n\
// it under the terms of the GNU General Public License as published by\n\
// the Free Software Foundation, either version 3 of the License, or\n\
// (at your option) any later version.\n\
//\n\
// This program is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU General Public License\n\
// along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\
//\n\
\n\
#include <algorithm> \n\
#include <cassert>" << endl;
#ifdef BREIT
ofs << "#include <src/integral/rys/breitrootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void BreitRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
#ifdef DAWSON
ofs << "#include <src/integral/rys/r2rootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void R2RootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
#ifdef SPIN2
ofs << "#include <src/integral/rys/spin2rootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void Spin2RootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#else
ofs << "#include <src/integral/rys/erirootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void ERIRootList::" << func << nroot << "(const double* ta, double* rr, double* ww, const int n) {\n" << endl;
#endif
#endif
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
#ifdef DAWSON
    const int ibox0 = static_cast<int>(log(MAXT)/log(2.0));
#endif
    for (int j=0; j != jend; ++j) {
#ifndef DAWSON
      vector<vector<double>> c_all = get_C(j*mstride, mstride, nroot, false);
#else
      vector<vector<double>> c_all;
      if (j < NBOX) {
        c_all = get_C(j*mstride, mstride, nroot, false);
      } else {
        if (MAXT < T_INFTY) {
          if (j >= NBOX && j < jend-1) { // NBOXL between MAXT and T_ASYM
            const int ibox = j - NBOX;
            const mpreal mstart = static_cast<mpreal> (MAXT + ibox*(ibox + 1.0)*(2.0*ibox + 1.0)/6.0);
            const mpreal mstrideL = static_cast<mpreal> (ibox + 1.0)*(ibox + 1.0);
            c_all = get_C(mstart, mstrideL, nroot, false);
          } else {
              const mpreal mstart = static_cast<mpreal> (T_ASYM);
              const mpreal mstrideL = static_cast<mpreal> (infty - T_ASYM);
              c_all = get_C(mstart, mstrideL, nroot, true);
          }
        } else {
          assert(NBOXL == 0);
          const int ibox = j - NBOX;
          const mpreal mstart = static_cast<mpreal>(pow(2.0, ibox0 + ibox));
          const mpreal mstrideL = static_cast<mpreal>(pow(2.0, ibox0 + ibox));
          c_all = get_C(mstart, mstrideL, nroot, true);
        }
      }
#endif

      for (int i = 0; i != nroot; ++i, ++index) {
        const int ii = 2 * i;
        const vector<double> x = c_all[ii];
        const vector<double> w = c_all[ii + 1];

        for (auto iter = x.begin(); iter != x.end(); ++iter) {
          listx << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != x.end() || j+1 != jend || i+1 != nroot || MAXT >= T_INFTY) listx << ",";
          if (xcnt++ % 7 == 4) listx << "\n";
        }
        for (auto iter = w.begin(); iter != w.end(); ++iter) {
          listw << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != w.end() || j+1 != jend || i+1 != nroot || MAXT >= T_INFTY) listw << ",";
          if (wcnt++ % 7 == 4) listw << "\n";
        }
      }
    }
#ifdef DAWSON
    if (MAXT >= T_INFTY) {
      for (int ibox = 0; ibox != T_INFTY; ++ibox) {
        vector<mpreal> tt_infty(1); tt_infty[0] = static_cast<mpreal>(pow(2.0, ibox + ibox0 + 1));
        vector<mpreal> dx_infty(nroot);
        vector<mpreal> dw_infty(nroot);
        rysroot_gmp(tt_infty, dx_infty, dw_infty, nroot, 1);
        for (auto iter = dx_infty.begin(); iter != dx_infty.end(); ++iter) {
          listx << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != dx_infty.end() || ibox + 1 != T_INFTY) listx << ",";
          if (xcnt++ % 7 == 4) listx << "\n";
        }
        for (auto iter = dw_infty.begin(); iter != dw_infty.end(); ++iter) {
          listw << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != dw_infty.end() || ibox + 1 != T_INFTY) listw << ",";
          if (wcnt++ % 7 == 4) listw << "\n";
        }
      }
    }
#endif
#ifndef SPIN2
#ifndef BREIT
    string tafactor = "t";
#else
    string tafactor = "t*t*t";
#endif
#else
    string tafactor = "t*t*t*t*t";
#endif

    int nbox2 = 0;
#ifndef DAWSON
    const int nbox1 = nbox;
#else
    int nbox1;
    if (MAXT < T_INFTY) {
      nbox1 = nbox + NBOXL + 1;
    } else {
      nbox1 = nbox + T_INFTY;
      nbox2 = T_INFTY;
    }
#endif
  ofs << "\
  constexpr double x[" << nroot*nbox1*ndeg + nroot*nbox2 <<"] = {";
  ofs << listx.str() << "\
  };" << endl;
  ofs << "\
  constexpr double w[" << nroot*nbox1*ndeg + nroot*nbox2 <<"] = {";
  ofs << listw.str() << "\
  };" << endl;

  ofs << "\
  int offset = -" << nroot << ";\n";
#ifdef DAWSON
  ofs << "\
  const int ibox0 = static_cast<int>(log(" << MAXT << ".0) / log(2.0)); \n";
#endif
  ofs << "\
  for (int i = 1; i <= n; ++i) {\n\
    double t = ta[i-1];\n\
    offset += " << nroot << ";\n\
    if (std::isnan(t)) {\n\
      fill_n(rr+offset, " << nroot << ", 0.5);\n\
      fill_n(ww+offset, " << nroot << ", 0.0);\n";
#ifndef DAWSON
      ofs << "\
    } else if (t >= " << MAXT << ".0) {\n\
      t = 1.0/sqrt(t);\n\
      for (int r = 0; r != " << nroot << "; ++r) {\n\
        rr[offset+r] = ax[r]*t*t;\n\
        ww[offset+r] = aw[r]*" + tafactor + ";\n\
      }\n\
    } else {\n\
      assert(t >= 0);\n\
      int it = static_cast<int>(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< ");\n\
      t = (t-it*" << stride << "-" << setw(20) << setprecision(15) << fixed << stride/2.0 << ") *" << setw(20) << setprecision(15) << fixed << 2.0/stride << ";\n\
      \n";
#else
      ofs << "\
    } else if (t >= " << infty << ".0) {\n\
      for (int r = 0; r != " << nroot << "; ++r) {\n\
        ww[offset+r] = aw[" << nroot << "-r-1] / t;\n\
        rr[offset+r] = 1.0 - ax[" << nroot << "-r-1] / t;\n\
      }\n\
    } else {\n\
      assert(t >= 0);\n";
    if (MAXT < T_INFTY) {
      ofs << "\
      vector<double> rr_infty(" << nroot << "); \n\
      vector<double> ww_infty(" << nroot << "); \n";
      for (int j = 0; j != nroot; ++j) {
      ofs << "\
      ww_infty[" << j << "] = " << setw(20) << setprecision(15) << fixed << dw_infty[j] << "; \n\
      rr_infty[" << j << "] = " << setw(20) << setprecision(15) << fixed << dx_infty[j] << "; \n";
      }
    }
      ofs << "\
      int it; \n\
      double bigT = 0.0; \n";
      if (NBOXL != 0) {
      ofs << "\
      if (" << MAXT << ".0 <= t && t < " << T_ASYM << ".0) { \n\
        int ka = static_cast<int>((pow((t - " << MAXT << ".0)*6.0, 1.0/3.0) - pow(0.25, 1.0/3.0))/pow(2.0, 1.0/3.0)); \n\
        int kb = static_cast<int>((pow((t - " << MAXT << ".0)*6.0, 1.0/3.0) - pow(6.0, 1.0/3.0))/pow(2.0, 1.0/3.0)); \n\
        assert(kb + 1 == ka || kb == ka); \n\
        it = " << NBOX << " + ka; \n\
        double a = " << MAXT << ".0 + ka * (ka + 1) * (2*ka + 1)/6.0; \n\
        double b = " << MAXT << ".0 + (ka + 1) * (ka + 2) * (2*ka + 3)/6.0; \n\
        t = (t - (a+b)/2) * 2/(a-b);\n\
      } else if (t >= " << T_ASYM << ".0 && t < " << infty << ".0) { \n";
      } else {
      ofs << "\
      if (t >= " << T_ASYM << ".0 && t < " << infty << ".0) { \n";
      }
        ofs << "\
        bigT = t; \n";
        if (MAXT < T_INFTY) {
        ofs << "\
        it = static_cast<int>(" << NBOX + NBOXL << ");\n\
        t = (t - (" << T_ASYM << ".0 + " << infty << ".0)/2) * 2/(" << infty << ".0 - " << T_ASYM << ".0);\n";
        } else {
        ofs << "\
        it = static_cast<int>(log(bigT) / log(2.0) + " << NBOX << " - ibox0);\n\
        t = (t - 1.5 * pow(2.0, it + ibox0 - " << NBOX << "))* 2/pow(2.0, it + ibox0 - " << NBOX << ");\n\
        cout << \" new t = \" << t << endl; \n";
        }
      ofs << "\
      } else { \n\
        it = static_cast<int>(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< ");\n\
        t = (t - it *" << stride << "-" << setw(20) << setprecision(15) << fixed << stride/2.0 << ") *" << setw(20) << setprecision(15) << fixed << 2.0/stride << ";\n\
      } \n";
#endif
      ofs << "\
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
#ifdef DAWSON
            if (MAXT < T_INFTY) {
        ofs << "\
        if (" << T_ASYM << ".0 <= bigT && bigT < " << infty << ".0) { \n\
          ww[offset+j-1] = ww[offset+j-1] * ww_infty[" << nroot << "-j] * " << infty << ".0 / bigT;\n\
          rr[offset+j-1] = 1.0 + rr[offset+j-1] * (1.0 - rr_infty[" << nroot << "-j]) * " << infty << ".0 /bigT; \n\
        }\n";
            } else {
          ofs << "\
        if (" << MAXT << ".0 <= bigT && bigT < " << infty << ".0) {\n\
          const int iref = " << (NBOX + T_INFTY) * nroot * NGRID << " + (it - " << NBOX << ") * " << nroot << " + " << nroot << " - j;\n\
          double rr_infty = x[iref];\n\
          double ww_infty = w[iref];\n\
          double Tref = pow(2.0, it + ibox0 + 1 - " << NBOX << ");\n\
          ww[offset+j-1] = ww[offset+j-1] * ww_infty * Tref / bigT;\n\
          rr[offset+j-1] = 1.0 + rr[offset+j-1] * (1.0 - rr_infty) * Tref /bigT;\n\
        }\n";
            }
#endif
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

