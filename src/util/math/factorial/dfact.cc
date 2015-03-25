//Compile with g++ -std=c++11 -L/opt/local/lib -lmpfr dfact.cc -o Dfact
#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpreal.h"
using namespace mpfr;
using namespace std;

mpreal dfac(const mpreal i) {
  if (i < 1.1) { return 1; }
  else { return i * dfac(i-2); }
}

int main() {
  mpreal::set_default_prec(10000);
  for (int i = 0; i != 18; ++i) {
    mpreal out = dfac(2*i-1);
    // out should be an integer
    if (abs(out-mpreal(out.toULong())) > 1.0e-10) assert(false);
    cout << "df_[" << setw(2) << i << "] = " << setw(20) << out.toULong() << "ull;";
    if (i % 4 == 3) cout << endl;
    else cout << " ";
  }
  cout << endl;
  return 0;
}
