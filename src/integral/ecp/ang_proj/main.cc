//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: January 31, 2014
//


#include "proj.h"

using namespace std;

int main() {


  cout << " Expansion of a Gaussian about a different centre " << endl;
  BesselI besselI;
  const int l = 10;
  const double x =0.1;
  cout << "I_n(x) where n = " << l << " and x = " << x << " is   " << besselI.compute(l, x) << endl;

  return 0;

}
