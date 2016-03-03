// Author:  Ryan D. Reynolds
// Date Created:  July 2013
// Based upon a related file for 1-D quadrature from April-May 2009

// These three lines all relate to the error function evaluation for complex numbers
constexpr double CUTOFF  = 0.5;    // If the argument's real part is smaller than this, use the Taylor expansion around 0; otherwise use the continued fraction
constexpr int NTAYLOR = 2000;    // Number of terms in the Taylor expansion used to represent the error function - 200 seems to be sufficient
constexpr int NTERMS  = 140000;  // Number of terms in the continued fraction used to represent the error function - 14,000 seems to be sufficient
constexpr int NBOYS   = 5000;    // Number of terms in the series expansion used for the Boys function at low values of T - 500 seems to be sufficient

// These next five lines all relate to the random number generation algorithm
constexpr int NTESTS  = 5000;  // This is the number of random T values we want to generate
constexpr double MAXERROR = 1.0e-15;  // If the error is greater than this, it will yell at us.
constexpr int MINTR =  -500;  // Minimum random T value - real component
constexpr int MAXTR =  1200;  // Maximum random T value - real component
constexpr int MINTI = -1000;  // Minimum random T value - imag component
constexpr int MAXTI =  1000;  // Maximum random T value - imag component

#include <fstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "mpreal.h"
#include "gmp_macros.h"
#include <algorithm>
#include <vector>
#include <cstdlib>                 // For debugging - includes a random number generator
#include <ctime>                   // For debugging - use time to generate a seed for the random number generator
#include <boost/lexical_cast.hpp>  // For debugging - used to interpret inputs to main{argc, **argv)

using namespace std;
using namespace mpfr;

// This function calculates the error function of complex number z using a continued fraction approximation.
// It is only valid when the real part of z is positive and larger than about 0.5 or so (depending on the number of terms used).
complex<mpreal> erffrac(const complex<mpreal> z){
  mpfr::mpreal::set_default_prec(GMPPREC);

  const mpreal one = "1.0";
  const mpreal half = "0.5";
  const mpreal sqrtpi = GMPPISQRT;
  complex<mpreal> K = one;
  complex<mpreal> result;
//  cout << "erffrac is used, with " << NTERMS << " terms." << endl;

  for(int j = NTERMS; j!=0; j--)
    K = half * j / (z + K);
  return one - exp(- z*z) / sqrtpi / (z + K) ;
}

// This function calculates the error function of complex number z using a Taylor series expansion centered at zero.
// It is only valid when z is close to zero.
complex<mpreal> erfsum(const complex<mpreal> z){
  mpfr::mpreal::set_default_prec(GMPPREC);

  const mpreal two = "2.0";
  const mpreal one = "1.0";
  const mpreal sqrtpi = GMPPISQRT;
  complex<mpreal> factor(0,0);
  complex<mpreal> term;

//  cout << "erfsum is used, with " << NTAYLOR << " terms." << endl;
  for(int j = 0; j != NTAYLOR; j++){
    mpreal power = two * j + one;
    mpreal fac = one;
    for (int k = j; k > 0; k--){
      fac *= k;
    }
    term = pow(-one,j) * pow(z,power) / fac / (two * j + one);
    factor += term;
  }
  return two / sqrtpi * factor;
}

// This is a generally-applicable function to calculate the error function of complex z.
// It calls whichever of the two above functions is more appropriate given the value of z.
complex<mpreal> complexerf(const complex<mpreal> z){
  mpfr::mpreal::set_default_prec(GMPPREC);

  const mpreal zero = "0.0";
  const mpreal cutoff = CUTOFF;

  const complex<mpreal> zz = z.real() < zero ? -z : z;

  complex<mpreal> erf;
  if (zz.real() < cutoff) {
    erf =  erfsum(zz);
  } else {
    erf =  erffrac(zz);
  }

  return z.real() < zero ? -erf : erf;
}

// Here is a function to calculate the zeroth order Boys function, given complex T.  It is only valid when the real part of
// T is negative and large.  See notebook page RDR-003-22 for details.
// It turns out that for problems of chemical interest, T should never be small enough to need this.
complex<mpreal> lowboys(const complex<mpreal> T){
  mpfr::mpreal::set_default_prec(GMPPREC);

  const mpreal two = "2.0";
  const mpreal one = "1.0";
  complex<mpreal> expansion = -one / (two * T);
  complex<mpreal> term;
  mpreal power;

  for(int k = 2; k != NBOYS + 2; k++){
    mpreal fac2 = one;
    int lim = 2*k - 3;
    for (int j = 1; j <= lim; j+=2){
      fac2 *= j;
    }
    power = k;
    term = pow((-one/(two*T)),power) * fac2;
    expansion += term;
  }
  return expansion * exp(-T);
}

// Here is a random number generator for the calculation of a random complex number
complex<mpreal> randcomplex(const int mintr, const int maxtr, const int minti, const int maxti) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  const int rrange = maxtr-mintr;
  const int irange = maxti-minti;

  const mpreal randreal = (rand()%(rrange*100));
  const mpreal randimag = (rand()%(irange*100));
  const mpreal realdouble = randreal/100 + mintr;
  const mpreal imagdouble = randimag/100 + minti;
  complex<mpreal> random (realdouble, imagdouble);
  return random;
}

// Here is a random number generator for the calculation of a random complex number, where the imaginary part cannot be greater than the real
complex<mpreal> randconstrained(const int mintr, const int maxtr) {
  mpfr::mpreal::set_default_prec(GMPPREC);

  const int rrange = maxtr-mintr;
  const int source = rand();
  const mpreal randreal = (source%(rrange*100));
  const mpreal realdouble = randreal/100 + mintr;
  const int randint = (source%(rrange*100));
  const int realint = (randint/100 + mintr);

  const int irange = 2*fabs(realint);
  const mpreal randimag = (rand()%(irange*100));
  const mpreal imagdouble = randimag/100 - fabs(realdouble);

  complex<mpreal> random (realdouble, imagdouble);
  return random;
}

void sortroots(vector<complex<mpreal>>& dx, vector<complex<mpreal>>& dw, const int nrank){
  mpfr::mpreal::set_default_prec(GMPPREC);
  bool sorted;
  do {
    sorted = true;
    for (int i = 0; i != nrank-1; i++) {
      if ( dx[i].real() > dx[i+1].real() ) {
        swap (dx[i],dx[i+1]);
        swap (dw[i],dw[i+1]);
        sorted = false;
      }
    }
  } while (sorted == false);
}

void sortweights(vector<complex<mpreal>>& dx, vector<complex<mpreal>>& dw, const int nrank){
  mpfr::mpreal::set_default_prec(GMPPREC);
  bool sorted;
  do {
    sorted = true;
    for (int i = 0; i != nrank-1; i++) {
      if ( fabs(dw[i]) < fabs(dw[i+1]) ) {
        swap (dx[i],dx[i+1]);
        swap (dw[i],dw[i+1]);
        sorted = false;
      }
    }
  } while (sorted == false);
}

// Here is the meat and potatoes of this file:  A function to determine the roots and weights for the
// Gaussian quadrature of an electron repulsion integral, given the value of T.
void complex_rysroot_gmp(const complex<mpreal>& ta, vector<complex<mpreal>>& dx, vector<complex<mpreal>>& dw, const int nrank){

  mpfr::mpreal::set_default_prec(GMPPREC);
  complex<mpreal> mlp[40];
  complex<mpreal> sigma[40];
  complex<mpreal> fm[40];
  complex<mpreal> w[40];
  complex<mpreal> x[40];
  complex<mpreal> lp[40];

  complex<mpreal> a[40];
  complex<mpreal> b[40];

  const complex<mpreal> T = ta;
  complex<mpreal> mone;

  {
    const mpreal zero = "0.0";
    const mpreal half = "0.5";
    const mpreal one = "1.0";
    const mpreal two = "2.0";

    const complex<mpreal> sqrtt = sqrt(T);

    const mpreal pi = GMPPI;
    const mpreal sqrtpi = GMPPISQRT;
    const complex<mpreal> halfpT = half / T;


    // First, we have to calculate the Boys function at zeroth order...
    const mpreal cutoff2 = "-500.0";
    if (T.real() > cutoff2){
      fm[0] = sqrt(pi) / sqrtt * half * complexerf(sqrtt);
    } else {
      fm[0] = lowboys(T);
    }

    // This Recurrence relation gives the solution to the Boys function of a given order
    // fm[i] = \int_0^1 e^{Tt}   / (2*sqrt(t)) t^i    dt    in Toru's preferred expression
    // fm[i] = \int_0^1 e^{Tt^2}               t^{2i} dt    using the definitions from Flocke's paper
    for (int i = 1; i != 40; ++i) {
      fm[i] = halfpT * ((two*i-one) * fm[i-1] - exp(-T));

      assert ( fabs(fm[i].real() < 1) || T.real() < 0 );
      // If this assertion has been triggered, than an unbelievable value for fm[i] was obtained.
      // This often results from inaccuracy in the complex error function that accumulates later on...  Try increasing the number of terms in the complexerf expansion.

    }


///*
//////////////////////////////////////////////////////////////
// To print off the values of the Boys function for a given T
// for(int i = 0; i!=2*(nrank+1); i++) cout << "\nF_" << i << "(T) = " << setprecision(15) << fm[i];
/////////////////////////////////////////////////////////////
//*/

    mone = fm[0];

  }
  {
    // Chebyshev algorithm; VERY unstable
    // Here we calculate the values of the recursion coefficients a_i and b_i
///*
    const mpreal mpone = "1.0";
    const int n = nrank;

    mlp[0] = mpone / fm[0];
    x[0] = fm[1] * mlp [0];

    // The loop of k runs nroot/2 times, rounded down (never for n == 1)
    for (int k = 0; k <= n - 2; k += 2) {
      for (int l = k; l <= 2 * n - k - 3; ++l) {
        sigma[l + 1] = fm[l + 2] - x[k] * fm[l + 1] - w[k] * sigma[l + 1] ;
      }
      // This section defines the odd-numbered entries of x, w, and mlp
      mlp[k + 1] = mpone / sigma[k + 1];
      x[k + 1] = - fm[k + 1] * mlp[k] + sigma[k + 2] * mlp[k + 1];
      w[k + 1] = sigma[k + 1] * mlp[k];

      if(k != n - 2) {
        for (int l = k + 1; l <= 2 * n - k - 4; ++l)
          fm[l + 1] = sigma[l + 2] - x[k + 1] * sigma[l + 1] - w[k + 1] * fm [l + 1];

        // This section defines the even-numbered entries of x, w, and mlp (except for the zeroth entry)
        mlp[k + 2] = mpone / fm[k + 2];
        x[k + 2] = -sigma[k + 2] * mlp[k + 1] + fm[k + 3] * mlp[k + 2];
        w[k + 2] = fm[k + 2] * mlp[k + 1];
      }
    }
//*/
////////////////////////////////////////////////////
/////  Setting up the Jacobi Matrix  ///////////////
////////////////////////////////////////////////////

// This portion of the code executes a four-term recurrence relation in order to find the elements of the Jacobi Matrix
// I think of the large phi array as an (n) by (2n) matrix.  Each row corresponds to a particular orthogonal polynomial p_i,
// and each column to a power of t.  The a_i and b_i values needed are simply the ratio of different phi elements.
// Each term of phi represents a weighted product of polynomials:  phi[2*n*i+k] = \int_0^1 e^-Tt/(2*sqrt{t}) t^k p_i p_i dt
/*
//    const int n = nrank;
    const int size = 2 * n * n;
    complex<mpreal> phi[size];

    // For the i = 0 case:
  {
    const mpreal zero = "0.0";
    const mpreal half = "0.5";
    const mpreal one = "1.0";

    const complex<mpreal> sqrtt = sqrt(T);

    const mpreal pi = GMPPI;
    const mpreal sqrtpi = GMPPISQRT;
    const complex<mpreal> halfpT = half / T;

    phi[0] = sqrt(pi) / sqrtt * half * complexerf(sqrtt);
        cout << "\n i = 0, k = 0 and the value " << setprecision(10) << phi[0];
    for (int i = 1; i != 2*n; ++i) {
      mpreal ii = i;
      phi[i] = halfpT * ((two*ii-one) * phi[i-1] - exp(-T));
        cout << "\n i = 0, k = " << i << " and the value " << setprecision(10) << phi[i];
    }
   }
//    for(int k=0; k!=2*n; k++){
//      phi[k] = fm[k];
//    }
      a[0] = phi[1] / phi[0];
    if (n > 1) {
      // For the i = 1 case:
      for(int k=0; k!=2*(n-1); k++){
        phi[2*n+k] = phi[k+2] - two * a[0] * phi[k+1] + a[0] * a[0] * phi[k];
        cout << "\n i = 1, k = " << k << " and the value " << setprecision(10) << phi[2*n+k] << " comes from " << phi[k+2] << phi[k+1] << phi[k];
      }
        a[1] = phi[2*n+1] / phi[2*n];
        b[1] = phi[2*n]   / phi[0];

      // For the i = 2 case:
        for(int k=0; k!=2*(n-2); k++){
          phi[4*n+k] = phi[2*n+k+2] - two * a[1] * phi[2*n+k+1] + a[1] * a[1] * phi[2*n+k] + b[1] * b[1] * phi[k];
          cout << "\n i = 2, k = " << k << " and the value " << setprecision(10) << phi[4*n+k] << " comes from " << phi[2*n+k+2] << phi[2*n+k+1]<< phi[2*n+k]<< phi[k];
        }
        a[2] = phi[4*n+1] / phi[4*n];
        b[2] = phi[4*n]   / phi[2*n];
////////////////////////////////////////////////////////
      // For all remaining values:  1 < i < n
      for(int i=3; i<n; i++){
        for(int k=0; k!=2*(n-i); k++){
          phi[2*n*i+k] = phi[2*n*(i-1)+k+2] - two * a[i-1] * phi[2*n*(i-1)+k+1] + a[i-1] * a[i-1] * phi[2*n*(i-1)+k] + b[i-1] * b[i-1] * phi[2*n*(i-2)+k];
          // The preceding line of code works correctly for i = 1.
          cout << "\n i = " << i << ", k = " << k << " and the value " << setprecision(10) << phi[2*n*i+k] << " comes from " << phi[2*n*(i-1)+k+2] << phi[2*n*(i-1)+k+1]<< phi[2*n*(i-1)+k]<< phi[2*n*(i-2)+k];
        }
        a[i] = phi[2*n*i+1] / phi[2*n*i];
        b[i] = phi[2*n*i]   / phi[2*n*(i-1)];
      }
    }
*////////////////////////////////////////////////////

    // Now we use what was done above to define the elements of dx and dw based on x and w.
    // dx is the same as x; we could probably simplify this as just dx = x
    // dw is a little more complicated - dw[i] is the square root of w[i+1]
    dx[0] = x[0];
    for (int i = 1; i !=n; ++i){
      dw[i - 1] = sqrt(w[i]);
      dx[i] = x[i];
    }
    dw[n - 1] = "0.0";
/*
    // For Ryan's alternate derivation of the a_i and b_i terms
    dx[0] = a[0];
    for (int i = 1; i !=n; ++i){
      dw[i - 1] = sqrt(b[i]);
      dx[i] = a[i];
    }
    dw[n - 1] = "0.0";
*/
//////////////////////////////////////////////////////////////
/*
// To print off the elements of the Jacobi matrix
cout << "\nToru's values: ";
for(int i = 0; i!=nrank; i++) cout << "\na_" << i << " = " << x[i];
for(int i = 1; i!=nrank; i++) cout << "\nb_" << i << " = " << w[i];
cout << "\nRyan's values: ";
for(int i = 0; i!=nrank; i++) cout << "\na_" << i << " = " << a[i];
for(int i = 1; i!=nrank; i++) cout << "\nb_" << i << " = " << b[i];
cout << "\nProducts used to generate elemetns: ";
for(int i = 0; i!=nrank; i++){
  cout << "\n";
  for(int j = 0; j!=2*nrank; j++){
    cout << setprecision(5) << phi[2*nrank*i+j] << ", ";
  }
}
//  cout << "\n";
//  cout << "\n";
//  for(int j = 0; j!=2*nrank*nrank; j++){
//    cout << setprecision(5) << phi[j] << ", ";
//  }
*/
/////////////////////////////////////////////////////////////
/*
    //////////////////////////////////////////
    ////  Printing out data for analysis  ////
    //////////////////////////////////////////
    cout << "Just before we get to 'solving the tri-diagonal linear equation'";
    cout << "\nThe elements of dx are: ";
    for(int i=0; i!=12; i++){
      cout << dx[i] << " ";
    }
    cout << "\nThe elements of dw are: ";
    for(int i=0; i!=12; i++){
      cout << dw[i] << " ";
    }
    cout << "\n";
    /////////////////////////////////////////
*/

///*
    // This next section finds the eigenvalues and eigenvectors of the tridiagonal matrix
    const mpreal zero = "0.0";
    const mpreal one = "1.0";
    const mpreal two = "2.0";

    lp[0] = one;
    for (int i = 1; i <= n + n - 2; ++i) lp[i] = zero; // Defining lp's elements as 1, 0, 0, 0, 0, 0...

    // Repeat this large for loop n times - Only one line of code remains after it
    // l identifies which root we're solving for, or something... each time, start at line1
    for (int l = 0; l <= n-1; ++l) {
      int iter = 0;

line1:
      int mm;
      for (mm = l; mm <= n - 2; ++mm) {
        mpreal dd = abs(dx[mm]) + abs(dx[mm + 1]);
        if (fabs(dw[mm]) + dd == dd) goto line2;
      }
      mm = n - 1;
line2:
      if(mm != l) {
        ++iter;
        if (iter == 100) throw logic_error("bad");
        complex<mpreal> g = (dx[l + 1] - dx[l]) / (dw[l] * two);
        complex<mpreal> rset;
        rset.real( sqrt(g.real() * g.real() + one) );
        rset.imag( g.imag() );
        const complex<mpreal> r = rset;
        // const complex<mpreal> r = sqrt(g * g + one);
        g = dx[mm] - dx[l] + dw[l] / ( g.real() >= zero ? ( g + r ) : ( g - r ) );
        complex<mpreal> s = one;
        complex<mpreal> c = one;
        complex<mpreal> p = zero;
        for (int i = mm - 1; i >= l; --i) {
          complex<mpreal> f = s * dw[i];
          complex<mpreal> bb = c * dw[i];
          // So f = bb at this point...
          complex<mpreal> r = sqrt(f * f + g * g);  // This appears to be distinct from the constant r declared above.
          dw[i + 1] = r;
          // If r = 0 (because f = g = 0), we execute the next three lines & are done with line2
          if(r == zero) {
            dx[i + 1] -= p; //Stylistic change here
            dx[mm] = zero;
            goto line1;
          }
          const complex<mpreal> pr = one / r;
          s = f * pr;
          c = g * pr;
          g = dx[i + 1] - p;
          complex<mpreal> td = c * bb;
          r = (dx[i] - g) * s + td * two;
          p = s * r;
          dx[i + 1] = g + p;
          g = c * r - bb;
          f = lp[i + 1];
          lp[i + 1] = s * lp[i] + c * f;
          lp[i] = c * lp[i] - s * f;
        }
        dx[l] -= p; //Stylistic change here
        dw[l] = g;
        dw[mm] = zero;
        goto line1;
      }
    }
//*/
    // We redefine dw at the end here...
    for (int i = 0; i != n; ++i) dw[i] = lp[i] * lp[i] * mone;
    for (int i = 0; i != n; i++) {
      if ( dx[i].real() < 0 || dx[i].real() > 1 || dx[i].imag() > 0) {
//        cout << "Unusual value for root " << i+1 << " of " << n << ";  T = " << T << ", root = " << dx[i] << ", weight = " << dw[i] << ".  " << endl;
//        assert ( fabs(T.imag()) > fabs(T.real()) || T.real() < 0 );
          // I had added an assertion here.  For the time being, I will deactivate it since I'm not sure that it's right.  For some absurd values of T I can trigger this
          // assertion without other things occuring that I know are wrong.
          // If this assertion has been triggered, then you obtained a root with an unexpected value, not between 0 and 1.  I have seen this
          // result from errors in the complex error function, but don't know of any examples that would cause problems here and not in the Boys function assert statement above.
      }
    }
  }
//if (nrank == 2) sortweights(dx, dw, nrank);
//if (nrank > 2) sortroots(dx, dw, nrank);
sortroots(dx, dw, nrank);
}


////////////////////////////////////////////////
//// To print out the data for playtesting  ////
////////////////////////////////////////////////

#if 0
int main(int argc, char** argv){
using namespace boost;
  mpfr::mpreal::set_default_prec(GMPPREC);

///*
  // This next section allows you to specify the form of the integral you want to evaluate at runtime
  // The first three arguments of int(main) are the number of roots to be calculated, then the real and imaginary parts of T
  if(argc > 2) {
    const string nvalue = argv[1];
    const int n = lexical_cast<int>(nvalue);

    const string Treal = argv[2];
    const double Tr = lexical_cast<double>(Treal);

    const complex<mpreal> zero (0,0);
    complex<mpreal> T = zero;
    T.real(Tr);

    vector<complex<mpreal>> dx(50, zero);
    vector<complex<mpreal>> dw(50, zero);

    cout << endl << "For the " << n << "-point quadrature of an ERI where T = " << Tr;
    if(argc > 3){
      const string Timag = argv[3];
      const double Ti = lexical_cast<double>(Timag);
      T.imag(Ti);
      cout << " + " << Ti << "i";
    }
    cout << endl << "and W(t) is given by e^(-Tt) / (2*sqrt(t))" << endl;

    // Run rysroot function to generate roots and weights
    if (n>0) complex_rysroot_gmp(T, dx, dw, n);

    // Print out roots and weights obtained
    cout << "The roots are:   " << endl;
    for(int i=0; i!=n; i++){
      cout << setprecision(15) << fixed << dx[i] << " ";
    }
    cout << endl << "The weights are: " << endl;
    for(int i=0; i!=n; i++){
      cout << setprecision(15) << fixed << dw[i] << " ";
    }
//    cout << endl;
    const complex<mpreal> expT = exp(-T);
    cout << endl << "The SCALED roots are: " << endl;
    for(int i=0; i!=n; i++){
      cout << setprecision(15) << fixed << dx[i]/expT << " ";
    }
    cout << endl;
    cout << endl << "The SCALED weights are: " << endl;
    for(int i=0; i!=n; i++){
      cout << setprecision(15) << fixed << dw[i]/expT << " ";
    }
    cout << endl;

    // If you enter more than 3 parameters, the rest will be interpreted as coefficients of a polynomial (real then imaginary parts)
    // It will then calculate the definite integral from 0 to 1 of e^-Tt/(2 sqrt(t)) * {Your polynomial of t} dt
    if(argc > 4 && argc % 2 == 0){

      // Define the coefficients of the polynomial
      complex<mpreal> coefficient[argc-4];
      cout << "Using the polynomial P(t) = ";
      for(int i = argc-6; i >=0; i=i-2){
        int power = i / 2;
        coefficient[power].real(argv[i + 4]);
        coefficient[power].imag(argv[i + 5]);
        cout << "(" << coefficient[power].real() << " + " << coefficient[power].imag() << "i)t^" << setprecision(1) << power;
        if (i != 0)
          cout << " + ";
        else
          cout << endl;
      }

      // Evaluate the integral by Gaussian quadrature
      complex<mpreal> intterm;
      complex<mpreal> integral (0,0);
      for(int j = 0; j!=n; j++){                 // For each root and weight pair...
        complex<mpreal> polyterm;
        complex<mpreal> polynomial(0,0);
        for(int i = argc-6; i >=0; i -= 2){       // For each term in the polynomial
          mpreal power = i/2;
          int coef = i/2;
          polyterm = coefficient[coef] * pow(dx[j],power);
          polynomial = polynomial + polyterm;
        }
        intterm = dw[j] * polynomial;
        integral = integral + intterm;
      }
      cout << "The integral from 0 to 1 of W(t) * P(t) is " << setprecision (15) << integral << endl;
    }
  }
//*/
  if(argc > 1) {
    const string toggle = argv[1];

    // Here is a feature to facilitate checking of the complex error function and total integral quadrature.  It creates many random complex numbers,
    // evaluates them, and prints off the results in a format that makes it easy to check against Wolfram-Alpha or some other external source.
    if (toggle == "-w") {
      const complex<mpreal> zero (0,0);
      constexpr int runs = NTESTS;
      const int npoly = 2;
      complex<mpreal> values[runs];
      complex<mpreal> roots[runs];
      complex<mpreal> results[runs];
      complex<mpreal> integral[runs];
      srand (time(NULL));

      for (int i = 0; i!=runs; i++){
        // For checking of the error function
        values[i] = randcomplex(MINTR, MAXTR, MINTI, MAXTI);   // For unconstrained imaginary parts
//        values[i] = randconstrained(MINTR, MAXTR);             // For imaginary parts no greater than the real parts (in terms of abs)
        roots[i] = sqrt(values[i]);
        results[i] = complexerf(roots[i]);
//        results[i] = complexerf(values[i]);

        int order[npoly];
        int rank = 0;
        complex<mpreal> coef[npoly];
        vector<complex<mpreal>> dx(50, zero);
        vector<complex<mpreal>> dw(50, zero);
        complex<mpreal> intterm;

        for (int k = 0; k != npoly; k++){
          order[k] = rand()%30;
          if (order[k] > rank) rank = order[k];
          coef[k] = randcomplex(MINTR, MAXTR, MINTI, MAXTI);
        }

        const int nroot = rank/2+1;
        complex_rysroot_gmp(values[i], dx, dw, nroot);
        integral[i] = zero;
        for(int j = 0; j!=nroot; j++){                 // For each root and weight pair...
          complex<mpreal> polyterm;
          complex<mpreal> polynomial(0,0);
          for(int k = 0; k != npoly; k++){       // For each term in the polynomial
            const mpreal power = order[k];
            polyterm = coef[k] * pow(dx[j],power);
            polynomial = polynomial + polyterm;
          }
          intterm = dw[j] * polynomial;
          integral[i] += intterm;
        }

        cout << setprecision(8) << "T = " << values[i].real() << " + " << values[i].imag() << " i,  sqrt(T) = " << roots[i].real() << " + " << roots[i].imag() << " i" << endl;
        cout << "rank = " << rank << "; nroot = " << nroot << endl;
//        cout << "erf( sqrt( " << values[i].real() << " + i * " << values[i].imag() << ") ) to 15 digits" << endl;
//////        cout << "erf( " << values[i].real() << " + i * " << values[i].imag() << ") to 15 digits" << endl;
//        cout << "= " << setprecision(15) << results[i].real() << " + " << results[i].imag() << " i" << endl;
        cout << "integral from 0 to 1 of  e^(-(" << values[i].real() << " + i * " << values[i].imag() << ")t)/(2*sqrt(t)) ( ";
        for (int k = 0; k != npoly; k++){
          cout <<  "(" << coef[k].real() << " + i * " << coef[k].imag() << ") * t^" << order[k] << " + ";
        }
        cout << "0 ) dt to 15 digits" << endl;
///*
        complex<mpreal> coef2[60];
        for (int q = 0; q != 30; q++){
          for (int k = 0; k != npoly; k++){
            if (q == order[k]) {
              coef2[q] = coef[k];
            }
          }
        }
        cout << "./gen " << nroot << " " << values[i].real() << " " << values[i].imag() << " ";
        for (int k = 0; k<=rank; k++){
          cout << coef2[k].real() << " " << coef2[k].imag() << " ";
        }
        cout << endl;
//*/
        cout << "= " << setprecision (15) << integral[i].real() << " + " << integral[i].imag() << "i" << endl << endl;
      }

/*
      ofstream ofs;
      string filename = "erftestvalues.txt";
      ofs.open(filename.c_str());
      for (int i = 0; i!=runs; i++){
        ofs << "erf( " << values[i].real() << " + i * " << values[i].imag() << " ) to 15 digits" << endl;
        ofs << values[i].real() << " , " << values[i].imag() << " , ";
        ofs << setprecision(15) << results[i].real() << " , " << results[i].imag() << " ,       ,       "<< endl << endl;
      }
    ofs.close();
*/
    }  // End -w test

    // Here is a general automated test:  Generate NTESTS random values of T, and for each one evalute the Boys function
    // of a random order both through the known expressions and by quadrature of the integral; compare results to ensure accuracy.
    // If we have a problem calculating roots or weights, this should reveal it.
    if (toggle == "-t") {
///*
      int failcount = 0;
      constexpr int runs = NTESTS;
      const mpreal half = "0.5";
      const mpreal one = "1.0";
      const mpreal two = "2.0";
      const mpreal pi = GMPPI;
      const mpreal sqrtpi = GMPPISQRT;
      const complex<mpreal> czero (0,0);

      srand (time(NULL));
      complex<mpreal> Trandom[runs];
      int order[runs];

      for (int i = 0; i != runs; i++){
        Trandom[i] = randcomplex(MINTR, MAXTR, MINTI, MAXTI);
        order[i] = (rand()%30);
        const int nroot = order[i]/2+1;
        cout << "T = " << Trandom[i] << ", order = " << order[i] << "\n";

        vector<complex<mpreal>> dx(50, czero);
        vector<complex<mpreal>> dw(50, czero);
        complex_rysroot_gmp(Trandom[i], dx, dw, nroot);

        // Evaluate the integral by Gaussian quadrature
        complex<mpreal> integral (0,0);
        mpreal power = order[i];
        for(int j = 0; j != nroot; j++){
          const complex<mpreal> polynomial = one * pow(dx[j],power);
          integral += dw[j] * polynomial;
        }
//        cout << "The integral from 0 to 1 of W(t) * P(t) is " << setprecision (15) << integral << "\n";

        const complex<mpreal> sqrtt = sqrt(Trandom[i]);
        const complex<mpreal> halfpT = half / Trandom[i];
        complex<mpreal> F_[40];
        F_[0] = sqrt(pi) / sqrtt * half * complexerf(sqrtt);
        for (int k = 1; k != 40; ++k) {
          F_[k] = halfpT * ((two*k - one) * F_[k-1] - exp(-Trandom[i]));
        }
//        cout << "The Boys function of order " << order[i] << " is " << F_[order[i]] << "\n";
        mpreal realerror = integral.real() - F_[order[i]].real();
        mpreal imagerror = integral.imag() - F_[order[i]].imag();
        mpreal realrelative = fabs(realerror / F_[order[i]].real());
        mpreal imagrelative = fabs(realerror / F_[order[i]].imag());
        mpreal maxerror = MAXERROR;
//        cout << "Error = " << error << "\n";
        if (fabs(realerror) > maxerror || fabs(imagerror) > maxerror) {
          cout << "Absolute error too high! T = " << Trandom[i] << " and nroot = " << order[i]/2+1 << "." << endl;
          cout << setprecision(15) << "Integral = " << integral << endl;
          cout << setprecision(15) << "Boys =     " << F_[order[i]] << endl;
          ++failcount;
        } else {
          if (realrelative > maxerror || imagrelative > maxerror) {
            cout << "Relative error too high! T = " << Trandom[i] << " and nroot = " << order[i]/2+1 << "." << endl;
            cout << setprecision(15) << "Integral = " << integral << endl;
            cout << setprecision(15) << "Boys =     " << F_[order[i]] << endl;
            ++failcount;
          } else {
            cout << "success" << endl;
          }
        }
      }     // End internal test for a given random number (for loop)
      cout << "Total fails: " << failcount << ", failure rate = " << failcount/NTESTS*100 << "%" << endl;
    }       // End -t argument
  }         // End 'if argc > 1'
//*/
  return 0;
}           // End main function
#endif
