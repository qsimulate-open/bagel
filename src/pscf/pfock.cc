//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <src/pscf/pfock.h>
#include <src/util/f77.h>
#include <cstring>

typedef std::shared_ptr<PFock> RefPFock;
typedef std::shared_ptr<PHcore> RefPHcore;
typedef std::shared_ptr<PGeometry> RefPGeometry;
typedef std::shared_ptr<PMatrix1e> RefPMatrix1e;

using namespace std;

PFock::PFock(const RefPGeometry g, const RefPFock prev, const RefPMatrix1e den, const vector<double>& shw, const int th, const bool dir)
 : PMatrix1e(g), previous_(prev), density_(den), schwarz_(shw), S2_(th), direct_(dir) {

  const int unit = 1;
  const complex<double> one(1.0, 0.0);

  // First construct the 2-e Fock matrices
  pfock_two_electron_part();

  PMatrix1e tmp = this->ft();
  zcopy_(&totalsize_, tmp.data()->front(), &unit, data_->front(), &unit); 

  // Finally, add the old Fock matrices
  zaxpy_(&totalsize_, &one, prev->data()->front(), &unit, data_->front(), &unit); 


}

PFock::PFock(const RefPGeometry g, const RefPFock prev, const RefPMatrix1e den, const vector<double>& shw, const int th, const bool dir, shared_ptr<PCompFile<ERIBatch> > fl)
 : PMatrix1e(g), previous_(prev), density_(den), schwarz_(shw), S2_(th), direct_(dir) {
  assert(!direct());

  file_ = fl;

  const int unit = 1;
  const complex<double> one(1.0, 0.0);
  // First construct the 2-e Fock matrices
  pfock_two_electron_part();
  PMatrix1e tmp = this->ft();
  zcopy_(&totalsize_, tmp.data()->front(), &unit, data_->front(), &unit); 
  // Finally, add the old Fock matrices
  zaxpy_(&totalsize_, &one, prev->data()->front(), &unit, data_->front(), &unit); 
}

PFock::PFock(const RefPGeometry g, const RefPHcore hcore) : PMatrix1e(g) {

  PMatrix1e hcore_m = hcore->ft();
  const int unit = 1;
  zcopy_(&totalsize_, hcore_m.data()->front(), &unit, data_->front(), &unit);

}


PFock::~PFock() {

}


