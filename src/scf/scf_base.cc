//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/scf_base.h>
#include <src/rysint/eribatch.h>
#include <src/util/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <boost/lexical_cast.hpp>

using namespace std;

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Matrix1e> RefMatrix1e;
typedef std::shared_ptr<Coeff> RefCoeff;
typedef std::shared_ptr<TildeX> RefTildeX;

SCF_base::SCF_base(std::multimap<std::string, std::string>& idat, const RefGeometry geom)
 : idata_(idat), geom_(geom), overlap_(new Overlap(geom)), hcore_(new Hcore(geom)) {

  eig_ = new double[geom_->nbasis()];
  hcore_->symmetrize();

  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  max_iter_ = read_input<int>(idata_, "maxiter_scf", max_iter_);
  thresh_overlap_ = read_input<double>(idata_, "thresh_overlap", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh_scf", thresh_scf_);
  string dd = read_input<string>(idata_, "diis", "gradient");
  if (dd == "gradient") {
    density_change_ = false;
  } else if (dd == "density") {
    density_change_ = true;
  } else {
    throw runtime_error("unrecongnized option for DIIS error vectors");
  }

  RefTildeX tildex_tmp(new TildeX(overlap_, thresh_overlap_));
  tildex_ = tildex_tmp;

  init_schwarz();
}


SCF_base::~SCF_base() {
  delete[] eig_;
}




void SCF_base::init_schwarz() {
  const vector<RefAtom> atoms = geom_->atoms(); 
  vector<RefShell> basis; 
  for (vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter) {
    const vector<RefShell> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
  }
  const int size = basis.size();

  schwarz_.resize(size * size);
  for (int i0 = 0; i0 != size; ++i0) {
    const RefShell b0 = basis[i0];
    for (int i1 = i0; i1 != size; ++i1) {
      const RefShell b1 = basis[i1];

      vector<RefShell> input;
      input.push_back(b1);
      input.push_back(b0);
      input.push_back(b1);
      input.push_back(b0);
      ERIBatch eribatch(input, 1.0);
      eribatch.compute();
      const double* eridata = eribatch.data();
      const int datasize = eribatch.data_size();
      double cmax = 0.0;
      for (int xi = 0; xi != datasize; ++xi, ++eridata) {
        const double absed = fabs(*eridata);
        if (absed > cmax) cmax = absed; 
      }
      schwarz_[i0 * size + i1] = cmax;
      schwarz_[i1 * size + i0] = cmax;
    }
  }
}
