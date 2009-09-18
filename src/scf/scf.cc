//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/scf.h>
#include <src/scf/scf_macros.h>
#include <src/rysint/eribatch.h>
#include <src/scf/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>

using namespace std;

typedef boost::shared_ptr<Geometry> RefGeometry;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Matrix1e> RefMatrix1e;
typedef boost::shared_ptr<Fock> RefFock;
typedef boost::shared_ptr<Coeff> RefCoeff;
typedef boost::shared_ptr<TildeX> RefTildeX;

SCF::SCF(const RefGeometry geom) : geom_(geom), overlap_(new Overlap(geom)), hcore_(new Hcore(geom)) {
  RefTildeX tildex_tmp(new TildeX(overlap_));
  tildex_ = tildex_tmp;

  eig_ = new double[geom_->nbasis()];
  hcore_->symmetrize();

  init_shwarz();
}


SCF::~SCF() {
  delete[] eig_;
}


void SCF::compute() {

  const string space3 = "   "; 
  const bool highest_level = geom_->level() == 0;
  string indent = "  ";
  for (int i = 0; i != geom_->level(); ++i) indent += "|";
  if (!highest_level) indent += "  ";

  RefFock previous_fock;
  {
    RefFock fock(new Fock(geom_, hcore_));
    previous_fock = fock;
   
    Matrix1e intermediate = *tildex_ % *fock * *tildex_;
    intermediate.diagonalize(eig_);
    Coeff new_coeff(*tildex_ * intermediate);
    RefMatrix1e guess_density(new Matrix1e(new_coeff.form_density_rhf()));

    aodensity_ = guess_density;
  }

  if (highest_level) {
    cout << indent << "=== Nuclear Repulsion===" << endl << indent << endl;
    cout << indent << fixed << setprecision(10) << setw(15) << geom_->nuclear_repulsion() << endl;
    cout << endl; 
  }
  cout << indent << "=== RHF iteration (" + geom_->basisfile() + ")===" << endl << indent << endl;

  // starting SCF iteration

  DIIS<Matrix1e> diis(5);
  RefMatrix1e densitychange = aodensity_; // assumes hcore guess...

  for (int iter = 0; iter != MAX_ITER_SCF; ++iter) {
    int start = ::clock();

    RefFock fock(new Fock(geom_, previous_fock, densitychange, shwarz_));
    previous_fock = fock;

    Matrix1e intermediate = *tildex_ % *fock * *tildex_;
    intermediate.diagonalize(eig_);
    RefCoeff new_coeff(new Coeff((*tildex_) * intermediate));
    coeff_ = new_coeff;
    RefMatrix1e new_density(new Matrix1e(coeff_->form_density_rhf()));

//  RefMatrix1e error_vector(new Matrix1e(*fock * *aodensity_ * *overlap_ - *overlap_ * *aodensity_ * *fock));
    RefMatrix1e error_vector(new Matrix1e(*new_density - *aodensity_));
    const double error = error_vector->rms();

    double energy = (*aodensity_ * *hcore_).trace() + geom_->nuclear_repulsion();
    for (int i = 0; i != geom_->nocc() / 2; ++i) energy += eig_[i];

    int end = ::clock();
    cout << indent << setw(5) << iter << setw(18) << fixed << setprecision(10) << energy << space3 
                                      << setw(15) << error << setw(15) << setprecision(2) << (end - start) * 1.0e-6 << endl; 

    if (error < SCF_THRESH) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    }


    RefMatrix1e diis_density = diis.extrapolate(make_pair(new_density, error_vector));

    RefMatrix1e dtmp(new Matrix1e(*diis_density - *aodensity_));
    densitychange = dtmp; 
    aodensity_ = diis_density;
  }
}


void SCF::print_eig() {
  for (int i = 0; i != geom_->nbasis(); ++i) 
//for (int i = 0; i != geom_->nocc() / 2; ++i) 
    cout << fixed << setw(15) << setprecision(10) << eig_[i] << " ";
  cout << endl;

}


void SCF::init_shwarz() {
  const vector<RefAtom> atoms = geom_->atoms(); 
  vector<RefShell> basis; 
  for (vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter) {
    const vector<RefShell> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
  }
  const int size = basis.size();

  shwarz_.resize(size * size);
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
      shwarz_[i0 * size + i1] = cmax;
      shwarz_[i1 * size + i0] = cmax;
    }
  }
}
