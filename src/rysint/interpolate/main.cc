//
// Author : Toru Shiozaki
// Date   : May 2009
//

#define NGRID 12
#define MAXT 64
#define SQRTPI2 0.886226925452758013649083741671
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "mpreal.h"
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
  for (int k = 0; k != n; ++k) {
     const mpreal y = mpfr::cos(GMPPI * static_cast<mpreal>(k + 0.5) / n); 
     out[k] = y;
  }
  return out;
}

vector<vector<double> > get_C(const double tbase, const double stride, int rank) {
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
    Tpoints[i] = stride*0.5*cheb[i] + Tp; 
  }

  vector<pair<vector<mpreal>, vector<mpreal> > > table_reserve;
  vector<mpreal> dx(rank);
  vector<mpreal> dw(rank);
  vector<mpreal> ttt(1);
  for (vector<mpreal>::const_iterator titer = Tpoints.begin(); titer != Tpoints.end(); ++titer) { 
    ttt[0] = *titer;
    rysroot_gmp(ttt, dx, dw, rank, 1);
    table_reserve.push_back(make_pair(dx, dw));
  }

  vector<vector<double> > c;
  for (int ii = 0; ii != rank; ++ii) {
    vector<double> tc(n);
    vector<double> tc2(n);

    vector<mpreal> cdx, cdw;
    for (int j = 0; j != n; ++j) {
      cdx.push_back((table_reserve[j].first)[ii]);
      cdw.push_back((table_reserve[j].second)[ii]);
    }  

    const mpreal two = "2.0";
    const mpreal fac = two / n; 
    const mpreal pi = GMPPI;
    for (int j = 0; j != n; ++j) {
      mpreal sum = 0.0;
      mpreal sum2 = 0.0;
      for (int k = 0; k != n; ++k) {
        sum += cdx[k] * cos(pi * j * (k + 0.5) / n);
        sum2 += cdw[k] * cos(pi * j * (k + 0.5) / n);
      }
      tc[j] = sum * fac;
      tc2[j] = sum2 * fac;
    }
    c.push_back(tc);
    c.push_back(tc2);
  }
  return c;
}

int main() {

  mpfr::mpreal pi = GMPPI;
  mpfr::mpreal::set_default_prec(GMPPREC);
  vector<double> nbox_(14);
  for (int nroot=1; nroot!=14; ++nroot) { 
    nbox_[nroot] = 32;
  } 

  for (int nroot=1; nroot!=14; ++nroot) { // this is the outer most loop.
    vector<double> aroot;
    vector<double> aweight;
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
      assert(fabs(dw[j]*sqrt(t) - dw[j+nroot]*sqrt(s)) < 1.0e-16);
      aroot.push_back((double)dx[j]*t);
      aweight.push_back((double)dw[j]*sqrt(t));
    }
#endif

    const int ndeg = NGRID;
    const int nbox = nbox_[nroot];
    const int jend = nbox;
    const double stride = (double)MAXT/(double)nbox;
    ofstream ofs; 
#ifndef BREIT
    const string func = "eriroot";
#else
    const string func = "breitroot";
#endif
    string filename = "_" + func + "_" + lexical_cast<string>(nroot) + ".f";
    ofs.open(filename.c_str());
    ofs << "\
!/\n\
!/ Author : Toru Shiozaki\n\
!/ Machine generated code\n\
!/" << endl;
    ofs << "\
      subroutine " + func + lexical_cast<string>(nroot) + "(ta, rr, ww, n)\n\
      implicit none\n\
      integer i, j, n, offset, it, boxof\n\
      double precision t, t2, d, e, f, g\n\
      double precision ta(*), rr(*), ww(*)\n\
      double precision ax("<<nroot<<")\n\
      double precision aw("<<nroot<<")\n\
      data (ax(i), i = 1, "<<nroot<<")/\n";
    for (int j=0; j!= nroot; ++j) {
      ofs << "\
     &   " << scientific << setprecision(15) << setw(20) << aroot[j];
      if (j != nroot-1) ofs << ",";
      ofs << endl;
    }
    ofs << "\
     &/\n\
      data (aw(i), i = 1, "<<nroot<<")/\n";
    for (int j=0; j!= nroot; ++j) {
      ofs << "\
     &   " << scientific << setprecision(15) << setw(20) << aweight[j]; 
      if (j != nroot-1) ofs << ",";
      ofs << endl;
    }
    ofs << "\
     &/\n\
      double precision x("<<nroot*nbox*ndeg<<")\n\
      double precision w("<<nroot*nbox*ndeg<<")\n";

////////////////////////////////////////
// now creates data
////////////////////////////////////////
    stringstream listx, listw;
    string indent("     & ");
    int nblock = 0;
    int index = 0;
    double tiny = 1.0e-100;
    for (int j=0; j!= jend; ++j) {
      vector<vector<double> > c_all = get_C(j*stride,stride,nroot);

      for (int i = 0; i != nroot; ++i, ++index) {
        const int ii = 2 * i;
        const vector<double> x = c_all[ii]; 
        const vector<double> w = c_all[ii + 1]; 

        for (vector<double>::const_iterator iter = x.begin(); iter != x.end(); ++iter) {
          listx << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) << 
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != x.end() || j+1 != jend || i+1 != nroot) listx << ",";
          listx << "\n";
        }
        for (vector<double>::const_iterator iter = w.begin(); iter != w.end(); ++iter) {
          listw << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
          if (iter + 1 != w.end() || j+1 != jend || i+1 != nroot) listw << ",";
          listw << "\n";
        }
      }
    }
    ofs << "\
      data x /\n";
    ofs << listx.str() << "\
     &/" << endl;
    ofs << "\
      data w /\n";
    ofs << listw.str() << "\
     &/" << endl;

    ofs << "\
      offset = -" << nroot << "\n\
      do i=1, n\n\
        t = ta(i)\n\
        offset = offset + " << nroot << "\n\
        if (t < 0.0d0) then\n\
          rr(offset+1:offset+" << nroot << ") = 0.5d0\n\
          ww(offset+1:offset+" << nroot << ") = 0.0d0\n\
        else if (t >= " << MAXT << ".0d0) then\n\
          t = 1.0d0/dsqrt(t)\n\
          rr(offset+1:offset+" << nroot << ") = ax(1:" << nroot << ")*t*t\n\
          ww(offset+1:offset+" << nroot << ") = aw(1:" << nroot << ")*t\n\
        else\n\
          it = int(t*" << setw(20) << setprecision(15) << fixed << 1.0/stride<< "d0)\n\
          t = (t-it*" << stride << "-" << setw(20) << setprecision(15) << fixed << stride/2.0 << "d0)\n     &     *"
                      << setw(20) << setprecision(15) << fixed << 2.0/stride << "d0\n\
          t2 = t * 2.0d0\n\
          do j=1, " << nroot << "\n\
            boxof = it*" << ndeg*nroot << "+" << ndeg << "*(j-1)\n";
     assert((ndeg/2)*2 == ndeg);
     for (int i=ndeg; i!=0; --i) {
       if (i==ndeg) {
         ofs << "\
            d = x(boxof+" << i << ")\n\
            e = w(boxof+" << i << ")\n";
       } else if (i==ndeg-1) {
         ofs << "\
            f = t2*d + x(boxof+" << i << ")\n\
            g = t2*e + w(boxof+" << i << ")\n";
       } else if (i != 1 && ((i/2)*2 == i)) { // odd
         ofs << "\
            d = t2*f - d + x(boxof+" << i << ")\n\
            e = t2*g - e + w(boxof+" << i << ")\n";
       } else if (i != 1) { // even
         ofs << "\
            f = t2*d - f + x(boxof+" << i << ")\n\
            g = t2*e - g + w(boxof+" << i << ")\n";
       } else {
         ofs << "\
            rr(offset+j) = t*d - f + x(boxof+" << i << ")*0.5d0\n\
            ww(offset+j) = t*e - g + w(boxof+" << i << ")*0.5d0\n";
       }
     }
       
     ofs << "\
          enddo\n\
        endif\n\
      enddo\n\
      end subroutine" << endl;

    ofs.close();
  }


  return 0;
}

