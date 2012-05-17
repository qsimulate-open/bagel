//
// Author : Toru Shiozaki
// Date   : May 2012
//

#include <src/prop/dipole.h>
#include <iomanip>

using namespace std;

Dipole::Dipole(shared_ptr<const Matrix1e> d) : den_(d) {

}


Dipole::~Dipole() {

}


void Dipole::compute() {
  vector<double> out(3, 0.0);

  cout << "    * Permanent dipole moment: (" << setw(12) << setprecision(6) << out[0] << ", "
                                             << setw(12) << out[1] << ", " << setw(12) << out[2] << ")" << endl;
}
