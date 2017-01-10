//
// Author:  Ryan D. Reynolds
// Creation date:  July 29, 2013
// Interpolation file for Rys Quadrature - Complex arguments
//

//////////////////////////////////////////
//////Basic Setup and Namespaces//////////
//////////////////////////////////////////

//#define TESTING 20              // Define this to skip code generation and instead run the functions in the generated files.  The number you give defines the number of tests to be run.

constexpr double IMULT = 0.01;          // Used for the "-h" option, this defines the starting value of T.imag as a multiple of T.real
constexpr double MAXABS_ERROR = 1e-14;
constexpr double MAXREL_ERROR = 2e-14;

constexpr int RGRID1 = 14;         // Number of gridpoints to use for bins close to zero
constexpr int IGRID1 = 10;
constexpr int CUTOFF = 11;          // Number of bins to use a larger number of gridpoints for (along real axis)

constexpr int RGRID2 = 10;         // Number of gridpoints to use above CUTOFF
constexpr int IGRID2 = 8;

#include "../comperirootlist.h"  // Used only the for the -t testing/debugging architecture
#include <src/util/constants.h>
#include <iostream>      // Basic input/output for the keyboard/screen
#include <sstream>       // "String stream," allows you to use strings for input or something
#include <vector>        // For vectors, similar to arrays
#include <complex>       // For complex numbers!
#include <fstream>       // File stream, for writing to files
#include <cmath>         // A collection of pre-defined mathematical functions such as cos()
#include <string>        // Allows us to use strings as a variable type
#include <iomanip>       // Defines setw and setprecision, which I think just deal with how the output is presented
#include <map>           // Provides the "map" class templates.
#include "mpreal.h"      // Defines the mpreal data type and all its uses
#include "gmp_macros.h"  // This just defines the values of SETPREC and GMPPI
#include <boost/lexical_cast.hpp>  // This gives you the ability to convert between data types
#include <cassert>       // This provides a debugging tool so "assert" functions can be checked.  It shouldn't be needed in the final form?

using namespace std;
using namespace mpfr;
using namespace boost;
using namespace bagel;

//////////////////////////////////////////
//////Definitions of Functions////////////
//////////////////////////////////////////

// dsyev is defined externally, I think as part of the LAPACK package...?
// Its purpose is to obtain the eigenvalues and eigenvectors of a matrix in order to
//   calculate the aroot and aweight vectors; see below for details
extern "C" {
  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
}

// This function is also defined externally, in rys_gmp.cc.  See below for details.
void complex_rysroot_gmp(const complex<mpreal>& ta, vector<complex<mpreal>>& dx, vector<complex<mpreal>>& dw, const int nrank);

complex<double> randcomplexdoub(const double mintr, const double maxtr, const double minti, const double maxti) {
  const int rrange = (maxtr-mintr)*10000;      // Don't make the coefficient here larger than 10,000
  const int irange = (maxti-minti)*10000;

  const double randreal = (rand()%(rrange));
  const double randimag = (rand()%(irange));
  const double realdouble = randreal/10000 + mintr;
  const double imagdouble = randimag/10000 + minti;
  complex<double> random (realdouble, imagdouble);
  return random;
}


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

void complex_get_C (const complex<mpreal>& Tbase, const mpreal& Rstride, const mpreal& Istride, const int rank, const int RGRID, const int IGRID, vector<vector<double>>& cxr, vector<vector<double>>& cxi, vector<vector<double>>& cwr, vector<vector<double>>& cwi){
  mpfr::mpreal::set_default_prec(GMPPREC);
  using namespace std;
  cout << "TBase = " << Tbase << ", " << rank << " roots needed.  Using " << RGRID << " by " << IGRID << " gridpoints." << endl;
  const int nr = RGRID;
  const int ni = IGRID;
  const int ntot = nr * ni;

  ////////////////////////////////////////////////////////////////
  //// Generates Chebyshev nodes over the range Tmin to Tmax /////
  ////////////////////////////////////////////////////////////////

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const complex<mpreal> Tmin = Tbase;
  const complex<mpreal> Tmax ( Tbase.real() + Rstride , Tbase.imag() + Istride );
  const complex<mpreal> Tp = half * (Tmin + Tmax);  // Center of the Chebyshev interval
  vector<complex<mpreal>> Tpoints(ntot);
  vector<mpreal> chebR = chebft(nr);
  vector<mpreal> chebI = chebft(ni);
  for (int j = 0; j != ni; j++) {
    for (int i = 0; i != nr; i++) {
      Tpoints[j*nr+i].real( Rstride*half*chebR[i] + Tp.real() );
//      Tpoints[j*nr+i].imag( 0 );
      Tpoints[j*nr+i].imag( Istride*half*chebI[j] + Tp.imag() );
//      cout << "Generating Tpoints: # " << j*nr+i+1 << " = " << Tpoints[j*nr+i] << endl;
    }
  }

  vector<vector<vector<mpreal>>> table_reserve(ntot);
  for (int i = 0; i != ntot; i++){
    complex<mpreal> ttt = Tpoints[i];
    vector<complex<mpreal>> dx(rank);
    vector<complex<mpreal>> dw(rank);
    complex_rysroot_gmp(ttt, dx, dw, rank);
//    cout << "Testing Tpoint # " << i+1 << " = " << Tpoints[i] << endl;
    for (int j = 0; j!= rank; j++){
      vector<mpreal> vec(4);
      vec[0] = dx[j].real();
      vec[1] = dx[j].imag();
      vec[2] = dw[j].real();
      vec[3] = dw[j].imag();
      table_reserve[i].push_back(vec);
//      cout << "Pushing a column to table_reserve:  Gridpoint " << i+1 << ", Root " << j+1 << ", x = " << dx[j].real() << ", w = " << dw[j].real() << endl;
    }
  }

  for (int ii=0; ii != rank; ii++) {               // This for loop essentially reorganizes data.  We are taking the data from rysroot_gmp and
    vector<double> tcxr(ntot);                       // Ordering the 12 gridpoints for a given interpolation together, rather than having all
    vector<double> tcxi(ntot);                       // the roots & weights for a given Tpoint together
    vector<double> tcwr(ntot);
    vector<double> tcwi(ntot);

    vector<mpreal> cdxr, cdwr, cdxi, cdwi;

    for (int j = 0; j != ntot; j++) {
      cdxr.push_back(table_reserve[j][ii][0]);     // These contain the real & imag parts of the roots & weights at each Tpoint
      cdxi.push_back(table_reserve[j][ii][1]);     // They are ordered as such:
      cdwr.push_back(table_reserve[j][ii][2]);     // They run through all (144) Tpoints, giving the data for the first root/weight pair of each
      cdwi.push_back(table_reserve[j][ii][3]);     // Then we repeat the ntot=144-point cycle, reporting the second root/weight, and so on...
    }

    {
      // Now we check for convergence along the imaginary axis, and throw the data out.  We'll keep the data for the real axis.
      vector<double> tcxr2(ntot);
      vector<double> tcxi2(ntot);
      vector<double> tcwr2(ntot);
      vector<double> tcwi2(ntot);
      const mpreal two = "2.0";
      const mpreal half = "0.5";
      const mpreal fac = two / ni;
      const mpreal pi = GMPPI;
      for (int i = 0; i != nr; i++) {                 // For each gridpoint on real axis - noninteracting at this point
        for (int j = 0; j != ni; ++j) {               // For each tabulated data point on imaginary axis - Each loop creates value of tc__2
          mpreal sumxr = "0.0";
          mpreal sumxi = "0.0";
          mpreal sumwr = "0.0";
          mpreal sumwi = "0.0";
          for (int k = 0; k !=ni; ++k) {              // For each input root for interpolation along the imag axis - we're summing up the terms to get tc__
            sumxr += cdxr[k*nr + i] * mpfr::cos(pi * j * (k + half) / ni);
            sumxi += cdxi[k*nr + i] * mpfr::cos(pi * j * (k + half) / ni);
            sumwr += cdwr[k*nr + i] * mpfr::cos(pi * j * (k + half) / ni);
            sumwi += cdwi[k*nr + i] * mpfr::cos(pi * j * (k + half) / ni);
          }
          tcxr2[i*ni + j] = (sumxr * fac).toDouble() ;
          tcxi2[i*ni + j] = (sumxi * fac).toDouble() ;
          tcwr2[i*ni + j] = (sumwr * fac).toDouble() ;
          tcwi2[i*ni + j] = (sumwi * fac).toDouble() ;
        }

        assert (fabs(tcxr2[(i+1)*ni - 1]) < 1.0e-6 || fabs(tcwr2[(i+1)*ni - 1]) < 1.0e-6 || fabs(tcxi2[(i+1)*ni - 1]) < 1.0e-6 || fabs(tcwi2[(i+1)*ni - 1]) < 1.0e-6);
        // If this assertion has been triggered, then the polynomial interpolation along the imaginary axis gave a very poor fit for the values obtained for each gridpoint.
        // The fit was bad enough that it likely cannot be fixed by simply using more gridpoints.  Check the actual values obtained for each gridpoint; perhaps some
        // of them will be implausibly large or reported in a different order than for the other gridpoints.

        if (fabs(tcxr2[(i+1)*ni - 1]) > 1.0e-12 || fabs(tcwr2[(i+1)*ni - 1]) > 1.0e-12 || fabs(tcxi2[(i+1)*ni - 1]) > 1.0e-12 || fabs(tcwi2[(i+1)*ni - 1]) > 1.0e-12) {
        cout << "CHEB FAILED (imag axis) for root " << ii+1 << " of " << rank << ", T ranging from (" << Tpoints[i].real() << "," << Tbase.imag() << ") to (" << Tpoints[i].real() << "," << Tbase.imag() + Istride << ")" << endl;
          cout << setw(20) << Tpoints[(i+1)*ni-1] << "tcx/w2 " << setw(20) << tcxr2[(i+1)*ni - 1]  << setw(20) << tcxi2[(i+1)*ni-1]  << setw(20) << tcwr2[(i+1)*ni-1]  << setw(20) << tcwi2[(i+1)*ni-1]  << endl;
        } else {
//          cout << "cheb passed (imag axis) for root " << ii+1 << " of " << rank << ", T ranging from (" << Tpoints[i].real() << "," << Tbase.imag() << ") to (" << Tpoints[i].real() << "," << Tbase.imag() + Istride << ")" << endl;
        }
      }
    }

    // Creating two sums:  sum for cdx and sum2 for cdw
    // sum and sum2 are both calculated for each element j of tc and tc2
    // They add up the terms of cdx/w, each one multiplied by cos(pi*j*(k+1/2)/n)
    const mpreal two = "2.0";
    const mpreal half = "0.5";
    const mpreal fac = two / nr;
    const mpreal pi = GMPPI;
    for (int i = 0; i != ni; i++) {                 // For each gridpoint on imaginary axis - noninteracting at this point
      for (int j = 0; j != nr; ++j) {               // For each tabulated data point on real axis - Each loop creates value of tc__
        mpreal sumxr = "0.0";
        mpreal sumxi = "0.0";
        mpreal sumwr = "0.0";
        mpreal sumwi = "0.0";
        for (int k = 0; k !=nr; ++k) {              // For each input root for interpolation along the real axis - we're summing up the terms to get tc__
          sumxr += cdxr[i*nr + k] * mpfr::cos(pi * j * (k + half) / nr);
          sumxi += cdxi[i*nr + k] * mpfr::cos(pi * j * (k + half) / nr);
          sumwr += cdwr[i*nr + k] * mpfr::cos(pi * j * (k + half) / nr);
          sumwi += cdwi[i*nr + k] * mpfr::cos(pi * j * (k + half) / nr);
        }
        tcxr[i*nr + j] = (sumxr * fac).toDouble() ;
        tcxi[i*nr + j] = (sumxi * fac).toDouble() ;
        tcwr[i*nr + j] = (sumwr * fac).toDouble() ;
        tcwi[i*nr + j] = (sumwi * fac).toDouble() ;

/*
        cout << setprecision(25) << "Root   " << ii+1 << " for gridpoint " << i*nr+j+1 << ":" ;
        cout << "  cdx = (" << cdxr[j] << ", " << cdxi[j] << ")";
        cout << ", tcx = (" << tcxr[j] << ", " << tcxi[j] << ")" << endl;
        cout << "Weight " << ii+1 << " for gridpoint " << i*nr+j+1 << ":" ;
        cout << "  cdw = (" << cdwr[j] << ", " << cdwi[j] << ")";
        cout << ", tcw = (" << tcwr[j] << ", " << tcwi[j] << ")" << endl;
*/
      }
      // If the last value of tc__ is not tiny, output a warning and the data points?
      // This warning shows up if you use too few grid points in the interpolation!

      assert (fabs(tcxr[(i+1)*nr - 1]) < 1.0e-6 || fabs(tcwr[(i+1)*nr - 1]) < 1.0e-6 || fabs(tcxi[(i+1)*nr - 1]) < 1.0e-6 || fabs(tcwi[(i+1)*nr - 1]) < 1.0e-6);
      // If this assertion has been triggered, then the polynomial interpolation along the real axis gave a very poor fit for the values obtained for each gridpoint.
      // The fit was bad enough that it likely cannot be fixed by simply using more gridpoints.  Check the actual values obtained for each gridpoint; perhaps some
      // of them will be implausibly large or reported in a different order than for the other gridpoints.

      if (fabs(tcxr[(i+1)*nr - 1]) > 1.0e-12 || fabs(tcwr[(i+1)*nr - 1]) > 1.0e-12 || fabs(tcxi[(i+1)*nr - 1]) > 1.0e-12 || fabs(tcwi[(i+1)*nr - 1]) > 1.0e-12) {
        cout << "CHEB FAILED (real axis) for root " << ii+1 << " of " << rank << ", T ranging from (" << Tbase.real() << "," << Tpoints[i*nr].imag() << ") to (" << Tbase.real()+Rstride << "," << Tpoints[i*nr].imag() << ")" << endl;
          cout << setw(20) << Tpoints[(i+1)*nr-1] << "tcx/w " << setw(20) << tcxr[(i+1)*nr - 1]  << setw(20) << tcxi[(i+1)*nr-1]  << setw(20) << tcwr[(i+1)*nr-1]  << setw(20) << tcwi[(i+1)*nr-1]  << endl;
//        cout << " caution: cheb not converged  - Root " << ii+1 << setprecision(10) << fixed << " " << tcxr[(i+1)*nr - 1] << " "
//             << tcwr[(i+1)*nr - 1] << " " << tcxi[(i+1)*nr - 1] << " " << tcwi[(i+1)*nr - 1] << endl;
//        for (int p = 0; p !=nr; ++p) {
//          cout << setw(20) << Tpoints[i*nr + p] << "cdx/w " << setw(20) << cdxr[i*nr + p].toDouble() << setw(20) << cdxi[i*nr + p].toDouble()
//               << setw(20) << cdwr[i*nr + p].toDouble() << setw(20) << cdwi[i*nr + p].toDouble() << endl;
//        }
//        for (int p = 0; p !=nr; ++p) {
//          cout << setw(20) << Tpoints[i*nr + p] << "tcx/w " << setw(20) << tcxr[i*nr + p]  << setw(20) << tcxi[i*nr + p]  << setw(20) << tcwr[i*nr + p]  << setw(20) << tcwi[i*nr + p]  << endl;
//        }
      } else {
//        cout << "cheb passed (real axis) for root " << ii+1 << " of " << rank << ", T ranging from (" << Tbase.real() << "," << Tpoints[i*nr].imag() << ") to (" << Tbase.real()+Rstride << "," << Tpoints[i*nr].imag() << ")" << endl;
      }

    }
    cxr.push_back(tcxr);
    cxi.push_back(tcxi);
    cwr.push_back(tcwr);
    cwi.push_back(tcwi);
  }
}

//////////////////////////////////////////
//////  The "main" Function  /////////////
//////////////////////////////////////////

int main (int argc, char** argv) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  mpfr::mpreal pi = GMPPI;

#ifdef TESTING
  if(argc > 1) {
    const string toggle = argv[1];

    ///////////////////////////////////////////
    ///// Testing with random numbers -t  /////
    ///////////////////////////////////////////
    if (toggle == "-t") {

      // Here we define the number of T values we want to evaluate
      const int NT = TESTING;

      // Here we initialize the arrays that will contain our results
      const static ComplexERIRootList mapT;
      const complex<mpreal> czero (0,0);

      // Here we assign the values of each T we will evaluate
      complex<double> Ts[NT];
      const int seed = time(NULL);
      srand (seed);
      for (int i = 0; i != NT; i++) {
        const complex<mpreal> TTT = randcomplexdoub(-2, 120, -0.5, 0.5);    // Inputs are MINTR, MAXTR, MINTI, MAXTI in that order
        Ts[i].real( TTT.real().toDouble() );
        Ts[i].imag( TTT.imag().toDouble() );
      }
      double failcount[RYS_MAX];

      for (int nroot = 1; nroot <= RYS_MAX; nroot++){
        failcount[nroot-1] = 0;

        // Evaluate each one first using the mapping in the generated files
        complex<double> mapt[NT*nroot];
        complex<double> mapw[NT*nroot];
        complex<double> ryst[NT*nroot];
        complex<double> rysw[NT*nroot];
        mapT.root(nroot, Ts, mapt, mapw, NT);

        // Then by converting to mpreal and runing complex_rysroot
        for (int k = 0; k != NT; k++){
          const complex<mpreal> rysT = Ts[k];
          vector<complex<mpreal>> rystvec (50, czero);
          vector<complex<mpreal>> ryswvec (50, czero);
          complex_rysroot_gmp(rysT, rystvec, ryswvec, nroot);
          for (int j = 0; j != nroot; j++){
            ryst[k*nroot+j].real( rystvec[j].real().toDouble() );
            rysw[k*nroot+j].real( ryswvec[j].real().toDouble() );
            ryst[k*nroot+j].imag( rystvec[j].imag().toDouble() );
            rysw[k*nroot+j].imag( ryswvec[j].imag().toDouble() );
          }
        }

        // Automated accuracy check:
        for(int i = 0; i!=NT*nroot; i++){

          mpreal rootRerror = ryst[i].real() - mapt[i].real();
          mpreal rootIerror = ryst[i].imag() - mapt[i].imag();
          mpreal rootRrel = rootRerror / fabs(ryst[i]);
          mpreal rootIrel = rootIerror / fabs(ryst[i]);
          mpreal weightRerror = rysw[i].real() - mapw[i].real();
          mpreal weightIerror = rysw[i].imag() - mapw[i].imag();
          mpreal weightRrel = weightRerror / fabs(rysw[i]);
          mpreal weightIrel = weightIerror / fabs(rysw[i]);

          if (fabs(rootRerror) > MAXABS_ERROR || fabs(rootIerror) > MAXABS_ERROR || fabs(weightRerror) > MAXABS_ERROR || fabs(weightIerror) > MAXABS_ERROR) {
            cout << "Absolute error too high! T = " << Ts[i/nroot] << ", root " << i%nroot+1 << " of " << nroot << "." << endl;
            cout << setprecision(15) << "Mapped root  = " << mapt[i] << " and mapped weight  = " << mapw[i] << endl;
            cout << setprecision(15) << "Rysroot root = " << ryst[i] << " and rysroot weight = " << rysw[i] << endl;
            ++failcount[nroot-1];
          } else {
            if (rootRrel > MAXREL_ERROR || rootIrel > MAXREL_ERROR) {                                                                // To check relative accuracy of roots only
//            if (rootRrel > MAXREL_ERROR || rootIrel > MAXREL_ERROR || weightRrel > MAXREL_ERROR || weightIrel > MAXREL_ERROR) {    // To check relative accuracy of both roots and weights
              cout << "Relative error too high! T = " << Ts[i/nroot] << ", root " << i%nroot+1 << "of " << nroot << "." << endl;
              cout << setprecision(15) << "Mapped root  = " << mapt[i] << " and mapped weight  = " << mapw[i] << endl;
              cout << setprecision(15) << "Rysroot root = " << ryst[i] << " and rysroot weight = " << rysw[i] << endl;
              ++failcount[nroot-1];
//            } else {
//                cout << "success for T = " << Ts[i/nroot] << ", root " << i%nroot+1 << " of " << nroot << "." << endl;
            }
          }
        }
        cout << "Total fails for " << nroot << "-root case = " << failcount[nroot-1] << ", failure rate = " << setprecision(4) << failcount[nroot-1]/(NT*nroot)*100 << "%" << endl;
      }
      cout << "Random number seed = " << seed << ". " << endl;
      for (int nroot = 1; nroot <= RYS_MAX; nroot++) {
        cout << "Total fails for " << nroot << "-root case = " << failcount[nroot-1] << ", failure rate = " << setprecision(4) << failcount[nroot-1]/(NT*nroot)*100 << "%" << endl;
      }
    }   // End if argv[1] = -t

    ////////////////////////////////////////////
    ///// Limits of the high-T expansion -h  ///
    ////////////////////////////////////////////
    if (toggle == "-h") {

      double limit[RYS_MAX];
      for (int nroot = 1; nroot <= RYS_MAX; nroot++){

        // Here we define the number of T values we want to evaluate
        bool pass = true;
        const int NT = 1;

        // Here we initialize the arrays that will contain our results
        const static ComplexERIRootList mapT;
        const complex<mpreal> czero (0,0);
        complex<double> mapt[NT*nroot];
        complex<double> mapw[NT*nroot];
        complex<double> ryst[NT*nroot];
        complex<double> rysw[NT*nroot];

        // Here we assign the values of each T we will evaluate
        complex<double> Ts[NT];
        Ts[0].real(150.0);
        Ts[0].imag( Ts[0].real() * IMULT );

        while (pass == true){

          for (int j=0; j != NT; j++) {
            Ts[j] = Ts[j] * 0.985;
          }

          // Evaluate each one first using the mapping in the generated files
          mapT.root(nroot, Ts, mapt, mapw, NT);

          // Then by converting to mpreal and runing complex_rysroot
          for (int k = 0; k != NT; k++){
            const complex<mpreal> rysT = Ts[k];
            vector<complex<mpreal>> rystvec (50, czero);
            vector<complex<mpreal>> ryswvec (50, czero);
            complex_rysroot_gmp(rysT, rystvec, ryswvec, nroot);
            for (int j = 0; j != nroot; j++){
              ryst[k*nroot+j].real( rystvec[j].real().toDouble() );
              rysw[k*nroot+j].real( ryswvec[j].real().toDouble() );
              ryst[k*nroot+j].imag( rystvec[j].imag().toDouble() );
              rysw[k*nroot+j].imag( ryswvec[j].imag().toDouble() );
            }
          }

          // Automated accuracy check:
          for(int i = 0; i!=NT*nroot; i++){

            double rootRerror = ryst[i].real() - mapt[i].real();
            double rootIerror = ryst[i].imag() - mapt[i].imag();
            double rootRrel = fabs(rootRerror / ryst[i].real());
            double rootIrel = fabs(rootIerror / ryst[i].imag());
            double weightRerror = rysw[i].real() - mapw[i].real();
            double weightIerror = rysw[i].imag() - mapw[i].imag();
            double weightRrel = fabs(weightRerror / rysw[i].real());
            double weightIrel = fabs(weightIerror / rysw[i].imag());

            if (fabs(rootRerror) > MAXABS_ERROR || fabs(rootIerror) > MAXABS_ERROR || fabs(weightRerror) > MAXABS_ERROR || fabs(weightIerror) > MAXABS_ERROR) {
              pass = false;
              limit[nroot-1] = Ts[0].real();
              cout << setprecision(15) << "Root " << i%nroot+1 << " Absolute error fails!  (This one checked first.)" << endl;
            } else {
              if (rootRrel > MAXREL_ERROR || rootIrel > MAXREL_ERROR || weightRrel > MAXREL_ERROR || weightIrel > MAXREL_ERROR) {
                pass = false;
                limit[nroot-1] = Ts[0].real();
                cout << setprecision(15) << "Root " << i%nroot+1 << " Relative error fails!" << endl;
//              } else {
//              cout << "success for T = " << Ts[i/nroot] << ", root " << i%nroot+1 << " of " << nroot << "." << endl;
              }
            }
          }
//          cout << "For " << nroot << " roots, and T.imag = T.real * " << IMULT << ", evaluated at T.real of " << Ts[0].real() << endl;
        } // End while loop
//        cout << "For " << nroot << " roots, and T.imag = T.real * " << IMULT << ", the high-T approximation broke down at T.real of " << limit[nroot-1] << endl;
      }   // End for(nroot) loop

      for (int nroot = 1; nroot <= RYS_MAX; nroot++){
        cout << "For " << nroot << " roots, and T.imag = T.real * " << IMULT << ", the high-T approximation broke down at T.real of " << limit[nroot-1] << endl;
      }

    }     // End if argv[1] = -h
    if(argc > 2) {
      const string nvalue = argv[1];
      const int n = lexical_cast<int>(nvalue);

      const string Treal = argv[2];
      const double Tr = lexical_cast<double>(Treal);

      const complex<mpreal> zero (0,0);
      complex<mpreal> T = zero;
      T.real(Tr);

//      cout << endl << "For the " << n << "-point quadrature of an ERI where T = " << Tr;
      if(argc > 3){
        const string Timag = argv[3];
        const double Ti = lexical_cast<double>(Timag);
        T.imag(Ti);
      }

      // Run get_C to see what values it tabulates for the specified bin
      vector<vector<complex<double>>> C;
      const mpreal stride = "2.0";
//      if (n>0) C = complex_get_C (T, stride, n);   BROKEN DEBUG/CHECKING CODE

    }     // End if argc > 2 (Specifying Tbase)
  }       // End if argc > 1

#endif
#ifndef TESTING
/*///// TESTING REALLY ////////////
const complex<mpreal> Tbase(8,0);
const mpreal Rstride = "2.0";
const mpreal Istride = "2.0";
const int rankt = 4;
vector<vector<double>> cxr, cxi, cwr, cwi;
complex_get_C (Tbase, Rstride, Istride, rankt, cxr, cxi, cwr, cwi);
*////// DONE TESTING //////////////

  ////////////////////////////////////////////////////////
  //// Each loop here creates a file for a given rank ////
  ////////////////////////////////////////////////////////

  for (int nroot=1; nroot <= RYS_MAX; ++nroot) {

    // Here we define the range across which we will use the interpolation
    const int MAXTRlist[13] = {34,42,48,54,60,66,70,76,80,86,92,96,100};
    const int MAXTR = MAXTRlist[nroot-1];
    const int MINTR = -2;
    const int MINTI = 0;
    const double MAXTI = 0.5;

    const mpreal mrstride = "2.0";
    const mpreal mistride = "0.5";
    const double drstride = mrstride.toDouble();
    const double distride = mistride.toDouble();

    const int rbox  = ((MAXTR-MINTR)/drstride);         // Number of bins we're dividing the whole T.real range into
    const int ibox  = ((MAXTI-MINTI)/distride);         // Number of bins we're dividing the whole T.imag range into
    const int tbox = rbox*ibox;
    cout << "rbox = " << rbox << ", ibox = " << ibox << endl;
    assert( ((MAXTR-MINTR)/drstride)/rbox == 1 );       // This assertion ensures that the interpolation range is an integer multiple of drstride.
    assert( ((MAXTI-MINTI)/distride)/ibox == 1 );       // This assertion ensures that the interpolation range is an integer multiple of distride.

    cout << "Generating data for " << nroot << " roots." << endl;
    vector<double> aroot;
    vector<double> aweight;
    const int n = nroot*2;
    const int ntot1 = RGRID1 * IGRID1;
    const int ntot2 = RGRID2 * IGRID2;
    const int nval1 = nroot * ntot1;
    const int nval2 = nroot * ntot2;

    ////////////////////////////////////////////////////////////
    //// Calculating ax and aw for the high-T approximation ////
    ////////////////////////////////////////////////////////////

    // These next 13 lines get somewhat ugly but are not interesting
    // They are just here to define the inputs for the dsyev function
    // Google "dsyev" to see what they each mean
    double a[10000] = {0.0};
    double b[100];
    double c[500];
    for ( int i=0; i!=n; ++i) {
      a[i+i*n] = 0.0;
      if (i > 0) {
        const double ia = static_cast<double>(i);
        a[(i-1)+i*n] = ::sqrt(ia*0.5);             // Here we are setting the values of a 2nroot by 2nroot matrix
        a[i+(i-1)*n] = ::sqrt(ia*0.5);             // They're technically stored as a 1-D array though
      }
    }
    int nn = n*5;
    int info = 0;


    // Though it looks messy, dsyev only depends upon the number of roots; everything else just sets up the problem
    // We feed it an (nroot by nroot) matrix a, and it gives us the eigenvalues and eigenvectors of this matrix
    // After calling this function, b contains the eigenvalues and a the eigenvectors (one element at a time)
    dsyev_("v", "U", &n, a, &n, b, c, &nn, &info);

    // The elements of aroot are given by the squares of the eigenvalues
    // The elements of aweight are given by the sqares of the first element of each eigenvalue, times sqrt(pi)
    // There are half as many roots/weights as eigenvalues/eigenvectors, because some of the eigenvalues are negative & give the same answer
    for (int j = 0; j != nroot; ++j) {
      aroot.push_back(b[nroot+j]*b[nroot+j]);
      aweight.push_back(a[(nroot+j)*(nroot*2)]*a[(nroot+j)*(nroot*2)]*(sqrt(pi)).toDouble());
    }

    ///////////////////////////////////////////////////////////
    //// Now we generate the files for each value of nroot ////
    ///////////////////////////////////////////////////////////
    const string func = "complex_eriroot";
    ofstream ofs;
    string filename = "_" + func + "_" + lexical_cast<string>(nroot) + ".cc";
    ofs.open(("../" + filename).c_str());
    ofs << "\
//\n\
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: " + filename + "\n\
// Copyright (C) 2013 Toru Shiozaki\n\
//\n\
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>\n\
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
#include <algorithm>\n\
#include <complex>\n\
#include \"comperirootlist.h\"\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void ComplexERIRootList::" << func << nroot << "(const complex<double>* ta, complex<double>* rr, complex<double>* ww, const int n) {\n";

// These next two sections print the ax and aw vectors we generated based on
// the output of the "dsyev" function.
ofs << "\n\
  static constexpr double ax["<<nroot<<"] = {";
    for (int j=0; j!=nroot; ++j) {
      ofs << scientific << setprecision(15) << setw(20) << aroot[j];
      if (j!=nroot-1) ofs << ",";
      if (j%7 == 4) ofs << endl << "    ";
    }
    ofs << "};" << endl;
    ofs << "\
  static constexpr double aw["<<nroot<<"] = {";
    for (int j=0; j!=nroot; ++j) {
      ofs << scientific << setprecision(15) << setw(20) << aweight[j];
      if (j!=nroot-1) ofs << ",";
      if (j%7 == 4) ofs << endl << "    ";
    }
    ofs << "};" << endl;

    // Now we generate long strings (listx and listw) and dump all the x and w values into them to facilitate printing.
    // This is a simplification:  It forces you to use the same RBOX (currently 32) for all generated files
    // Toru's version defined a vector rbox_, so that you could use a different value depending on the rank of the polynomial
    stringstream listxr, listxi, listwr, listwi;
    string indent(" ");
    int xrcnt = 0;
    int xicnt = 0;
    int wrcnt = 0;
    int wicnt = 0;
    const double tiny = 1.0e-100;

    for (int k=0; k != ibox; ++k) {               // For each row of boxes along the imaginary axis
      for (int j=0; j != rbox; ++j) {             // For each box along the real axis

        const int jk = k*rbox+j;
        vector<vector<double>> cxr, cxi, cwr, cwi;
        complex<mpreal> Tbase(MINTR + j*mrstride, MINTI + k*mistride);
//        cout << "Tbase = " << Tbase << ", mrstride = " << mrstride << ", mistride = " << mistride << ", nroot = " << nroot << endl;

        if (j < CUTOFF) complex_get_C (Tbase, mrstride, mistride, nroot, RGRID1, IGRID1, cxr, cxi, cwr, cwi);
        else complex_get_C (Tbase, mrstride, mistride, nroot, RGRID2, IGRID2, cxr, cxi, cwr, cwi);

        for (int i = 0; i != nroot; ++i) {
          const vector<double> xr = cxr[i];
          const vector<double> xi = cxi[i];
          const vector<double> wr = cwr[i];
          const vector<double> wi = cwi[i];

          for (auto iter = xr.begin(); iter != xr.end(); ++iter) {
            listxr << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? "  " : " ") << setw(20) <<
                                (fabs(*iter) < tiny ? 0.0 : *iter);
            if (iter + 1 != xr.end() || jk+1 != tbox || i+1 != nroot) listxr << ",";
            if (xrcnt++ % 7 == 4) listxr << "\n";
          }

          for (auto iter = xi.begin(); iter != xi.end(); ++iter) {
            listxi << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? "  " : " ") << setw(20) <<
                                (fabs(*iter) < tiny ? 0.0 : *iter);
            if (iter + 1 != xi.end() || jk+1 != tbox || i+1 != nroot) listxi << ",";
            if (xicnt++ % 7 == 4) listxi << "\n";
          }

          for (auto iter = wr.begin(); iter != wr.end(); ++iter) {
            listwr << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? "  " : " ") << setw(20) <<
                                (fabs(*iter) < tiny ? 0.0 : *iter);
            if (iter + 1 != wr.end() || jk+1 != tbox || i+1 != nroot) listwr << ",";
            if (wrcnt++ % 7 == 4) listwr << "\n";
          }

          for (auto iter = wi.begin(); iter != wi.end(); ++iter) {
            listwi << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? "  " : " ") << setw(20) <<
                                (fabs(*iter) < tiny ? 0.0 : *iter);
            if (iter + 1 != wi.end() || jk+1 != tbox || i+1 != nroot) listwi << ",";
            if (wicnt++ % 7 == 4) listwi << "\n";
          }
        }
      }
    }

  // This next little chunk prints all the x and w values from their strings into the generated files
  // Note that I have removed a line that is unneded due to exclusion of BREIT & SPIN2:  string tafactor = "t";
    ofs << "\
  static constexpr double xr[" << ibox*(rbox*nval2+CUTOFF*(nval1-nval2)) << "] = {";
    ofs << listxr.str() << "\
  };" << endl;
    ofs << "\
  static constexpr double xi[" << ibox*(rbox*nval2+CUTOFF*(nval1-nval2)) << "] = {";
    ofs << listxi.str() << "\
  };" << endl;
    ofs << "\
  static constexpr double wr[" << ibox*(rbox*nval2+CUTOFF*(nval1-nval2)) << "] = {";
    ofs << listwr.str() << "\
  };" << endl;
    ofs << "\
  static constexpr double wi[" << ibox*(rbox*nval2+CUTOFF*(nval1-nval2)) << "] = {";
    ofs << listwi.str() << "\
  };" << endl;
    ofs << "\
";

  // At this point, we have printed off all the calculated x and w values.
  // The rest of this file is about code generation - producing some math that will appear in the generated files but not be run by main.cc.
  ofs << "\
  int offset = -" << nroot << ";\n\
  for (int i = 1; i <= n; ++i) {\n\
    complex<double> t = ta[i-1];\n\
    offset += " << nroot << ";\n\
    if (std::isnan(t.real())) {\n\
      fill_n(rr+offset, " << nroot << ", 0.5);\n\
      fill_n(ww+offset, " << nroot << ", 0.0);\n\
    } else if (t.real() < " << MINTR << ") {\n\
      throw runtime_error (\"ERROR!  Invalid T value!  Real part is too small.  Consider regenerating interpolation files with a larger domain or reducing the magnetic field strength.\");\n\
    } else if (t.real() >= " << MAXTR << ") {\n";
//      cout << \"T = \" << t << \", need " << nroot << " roots:  Used high-T approximation\" << endl;\n
      ofs << "\
      t = 1.0/sqrt(t);\n\
      for (int r = 0; r != " << nroot << "; ++r) {\n\
        rr[offset+r] = ax[r]*t*t;\n\
        ww[offset+r] = aw[r]*t;\n\
      }\n\
    } else if ( fabs(t.imag()) > " << MAXTI << "){\n\
      throw runtime_error (\"ERROR!  Invalid T value!  Magnitude of imaginary part is too large.  Consider regenerating interpolation files with a larger domain or reducing the magnetic field strength.\");\n\
    } else {\n\
      const complex<double> torig = t;\n\
      if (torig.imag() < 0) t = conj(torig);\n\
      int itr = static_cast<int>((t.real()-(" << MINTR <<"))*" << setw(20) << setprecision(15) << fixed << 1.0/drstride << ") ;\n\
      int iti = static_cast<int>((t.imag()-(" << MINTI <<"))*" << setw(20) << setprecision(15) << fixed << 1.0/distride << ") ;\n\
      double tr = (t.real()-itr*" << drstride << "-" << setw(20) << setprecision(15) << fixed << drstride/2.0 << " - " << MINTR << ") *" << setw(20) << setprecision(15) << fixed << 2.0/drstride << ";\n\
      double ti = (t.imag()-iti*" << distride << "-" << setw(20) << setprecision(15) << fixed << distride/2.0 << " - " << MINTI << ") *" << setw(20) << setprecision(15) << fixed << 2.0/distride << ";\n\
      const double tr2 = tr * 2.0;\n\
      const double ti2 = ti * 2.0;\n";
      for (int cycle = 0; cycle != 2; cycle++) {
        int RGRID, IGRID;
        if (cycle == 0) {
          RGRID = RGRID1;
          IGRID = IGRID1;
          ofs << "\
      if (itr < " << CUTOFF << ") {\n";
        } else {
          RGRID = RGRID2;
          IGRID = IGRID2;
          ofs << "\
      } else {\n";
        }
//        cout << \"T = \" << torig << \", need " << nroot << " roots: Used complex interpolation, with " << RGRID << " by " << IGRID << " gridpoints.\" << endl;
        ofs << "\
        for (int j=1; j <=" << nroot << "; ++j) {\n\
          double xrval[" << IGRID << "];\n\
          double xival[" << IGRID << "];\n\
          double wrval[" << IGRID << "];\n\
          double wival[" << IGRID << "];\n\
          for (int k=1; k<= " << IGRID << "; k++){\n";
        if (cycle == 0) {
          ofs << "\
            const int boxof = iti*" << (rbox * nval2 + CUTOFF * (nval1-nval2) ) << " + itr*" << nval1 << " + (j-1)*" << ntot1 << " + (k-1)*" << RGRID << ";\n";
        } else {
          ofs << "\
            const int boxof = iti*" << (rbox * nval2 + CUTOFF * (nval1-nval2) ) << " + itr*" << nval2 << " + " << CUTOFF*(nval1-nval2) << " + (j-1)*" << ntot2 << " + (k-1)*" << RGRID << ";\n";
        }
         assert((RGRID/2)*2 == RGRID);     // This assertion generates an error message if RGRID is odd
         assert((IGRID/2)*2 == IGRID);     // This assertion generates an error message if IGRID is odd
         for (int i=RGRID; i!=0; --i) {    // This "for" loop runs RGRID = 12 times
           if (i==RGRID) {
             // First time through:  i = RGRID = 12
             ofs << "\
            double dr = xr[boxof+" << i-1 << "];\n\
            double di = xi[boxof+" << i-1 << "];\n\
            double er = wr[boxof+" << i-1 << "];\n\
            double ei = wi[boxof+" << i-1 << "];\n";
           } else if (i==RGRID-1) {
             // Second time through:  i = RGRID-1 = 11
             ofs << "\
            double fr = tr2*dr + xr[boxof+" << i-1 << "];\n\
            double fi = tr2*di + xi[boxof+" << i-1 << "];\n\
            double gr = tr2*er + wr[boxof+" << i-1 << "];\n\
            double gi = tr2*ei + wi[boxof+" << i-1 << "];\n";
           } else if (i !=1 && ((i/2)*2 == i)) {
             // This next cyles runs for odd-numbered iterations (even values of i, because of the way integers round)
             ofs << "\
            dr = tr2*fr - dr + xr[boxof+" << i-1 << "];\n\
            di = tr2*fi - di + xi[boxof+" << i-1 << "];\n\
            er = tr2*gr - er + wr[boxof+" << i-1 << "];\n\
            ei = tr2*gi - ei + wi[boxof+" << i-1 << "];\n";
           } else if (i !=1) {
             // Now for even-numbered iterations (odd values of i)
             ofs << "\
            fr = tr2*dr - fr + xr[boxof+" << i-1 << "];\n\
            fi = tr2*di - fi + xi[boxof+" << i-1 << "];\n\
            gr = tr2*er - gr + wr[boxof+" << i-1 << "];\n\
            gi = tr2*ei - gi + wi[boxof+" << i-1 << "];\n";
           } else {
             // Last cycle:  i = 1
             ofs << "\
            xrval[k-1] = tr*dr - fr + xr[boxof+" << i-1 << "]*0.5;\n\
            xival[k-1] = tr*di - fi + xi[boxof+" << i-1 << "]*0.5;\n\
            wrval[k-1] = tr*er - gr + wr[boxof+" << i-1 << "]*0.5;\n\
            wival[k-1] = tr*ei - gi + wi[boxof+" << i-1 << "]*0.5;\n";
           }
         }
       ofs << "\
          }\n\
          const double denom = " << IGRID << ";\n\
          const double fac = 2 / denom;\n\
          const double pi = 3.141592653589793238462;\n\
          double tcxr[" << IGRID << "];\n\
          double tcxi[" << IGRID << "];\n\
          double tcwr[" << IGRID << "];\n\
          double tcwi[" << IGRID << "];\n\
          for (int b = 0; b != " << IGRID << "; ++b) {\n\
            double sumxr = 0;\n\
            double sumxi = 0;\n\
            double sumwr = 0;\n\
            double sumwi = 0;\n\
            double fac2 = pi * b / " << IGRID << ";\n\
            for (int k = 0; k !=" << IGRID << "; ++k) {\n\
              double fac3 = cos (fac2 * (k + 0.5) );\n\
              sumxr += xrval[k] * fac3;\n\
              sumxi += xival[k] * fac3;\n\
              sumwr += wrval[k] * fac3;\n\
              sumwi += wival[k] * fac3;\n\
            }\n\
            tcxr[b] = (sumxr * fac);\n\
            tcxi[b] = (sumxi * fac);\n\
            tcwr[b] = (sumwr * fac);\n\
            tcwi[b] = (sumwi * fac);\n\
          }\n";
         for (int i=IGRID; i!=0; --i) {    // This "for" loop runs IGRID = 12 times
           if (i==IGRID) {
             // First time through:  i = IGRID = 12
             ofs << "\
          double dr = tcxr[" << i-1 << "];\n\
          double di = tcxi[" << i-1 << "];\n\
          double er = tcwr[" << i-1 << "];\n\
          double ei = tcwi[" << i-1 << "];\n";
           } else if (i==IGRID-1) {
             // Second time through:  i = IGRID-1 = 11
             ofs << "\
          double fr = ti2*dr + tcxr[" << i-1 << "];\n\
          double fi = ti2*di + tcxi[" << i-1 << "];\n\
          double gr = ti2*er + tcwr[" << i-1 << "];\n\
          double gi = ti2*ei + tcwi[" << i-1 << "];\n";
           } else if (i !=1 && ((i/2)*2 == i)) {
             // This next cyles runs for odd-numbered iterations (even values of i, because of the way integers round)
             ofs << "\
          dr = ti2*fr - dr + tcxr[" << i-1 << "];\n\
          di = ti2*fi - di + tcxi[" << i-1 << "];\n\
          er = ti2*gr - er + tcwr[" << i-1 << "];\n\
          ei = ti2*gi - ei + tcwi[" << i-1 << "];\n";
           } else if (i !=1) {
             // Now for even-numbered iterations (odd values of i)
             ofs << "\
          fr = ti2*dr - fr + tcxr[" << i-1 << "];\n\
          fi = ti2*di - fi + tcxi[" << i-1 << "];\n\
          gr = ti2*er - gr + tcwr[" << i-1 << "];\n\
          gi = ti2*ei - gi + tcwi[" << i-1 << "];\n";
           } else {
             // Last cycle:  i = 1
             ofs << "\
          rr[offset+j-1].real( ti*dr - fr + tcxr[" << i-1 << "]*0.5 ) ;\n\
          ww[offset+j-1].real( ti*er - gr + tcwr[" << i-1 << "]*0.5 ) ;\n\
          if (torig.imag() < 0) {\n\
            rr[offset+j-1].imag(-1*( ti*di - fi + tcxi[" << i-1 << "]*0.5 )) ;\n\
            ww[offset+j-1].imag(-1*( ti*ei - gi + tcwi[" << i-1 << "]*0.5 )) ;\n\
          } else {\n\
            rr[offset+j-1].imag( ti*di - fi + tcxi[" << i-1 << "]*0.5 ) ;\n\
            ww[offset+j-1].imag( ti*ei - gi + tcwi[" << i-1 << "]*0.5 ) ;\n\
          }\n\
        }\n";
          }
        }
      }
ofs << "\
      }\n\
    }\n\
  }\n\
}";
    ofs.close();
//#endif
  }
/////////////////////////////////////////////
/////  End of file creation/nrank loop  /////
/////////////////////////////////////////////
#endif
  return 0;
}
