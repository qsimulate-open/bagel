#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpreal.h"
using namespace mpfr;
using namespace std;

mpreal fac(const mpreal i) {
  if (i < 1.1) { return 1; }
  else { return i * fac(i-1); }
}

int main() {
  mpreal::set_default_prec(10000);
  for (int i = 0; i != 21; ++i) {
    mpreal out = fac(i);
    // out should be an integer
    if (abs(out-mpreal(out.toULong())) > 1.0e-10) assert(false);
    cout << "f_[" << setw(2) << i << "] = " << setw(20) << out.toULong() << "ull;";
    if (i % 4 == 3) cout << endl;
    else cout << " ";
  }
  return 0;
}
