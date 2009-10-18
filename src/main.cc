//
// Author: Toru Shiozaki
// Date  : April 2009
//
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <src/pscf/pgeometry.h>
#include <src/pscf/pscf.h>
#include <src/pscf/pscf_disk.h>
#include <src/pmp2/pmp2.h>
#include <src/scf/geometry.h>
#include <src/scf/scf.h>
#include <src/global.h>

using namespace std;

int main(int argc, char** argv) {

  try {
    print_header();

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    const int depth_basis = count_string(input, "Basis");

    const bool periodic = (bool)count_string(input, "Periodic");
    const bool domp2 = (bool)count_string(input, "MP2");

    typedef boost::shared_ptr<Geometry> RefGeom;
    typedef boost::shared_ptr<PGeometry> RefPGeom;
    typedef boost::shared_ptr<PSCF_DISK> RefPSCF_DISK;

    RefPSCF_DISK pscf;

    for (int i = depth_basis - 1; i != -1; --i) {
      if (!periodic) {
        RefGeom geom(new Geometry(input, i));
        SCF scf(geom);
        scf.compute();
      } else {
        RefPGeom pgeom(new PGeometry(input, i));
        RefPSCF_DISK tmp(new PSCF_DISK(pgeom));
        pscf = tmp;
        pscf->compute();
      }
    }

    if (periodic && domp2) {
      PMP2 pmp2(pscf->geom(), pscf->coeff(), pscf->eig(), pscf->ao_eri());
      pmp2.compute();
    }
  } catch (const exception &e) {
    cout << "Caught : " << e.what() << endl;
    cout << "Type   : " << typeid(e).name() << endl;
  }

  return 0;
}

