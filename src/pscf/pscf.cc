//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <src/pscf/pscf.h>
#include <src/pscf/pcoeff.h>
#include <src/scf/scf_macros.h>
#include <src/pscf/pdiis.h>
#include <src/rysint/eribatch.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<PGeometry> RefPGeometry;
typedef boost::shared_ptr<PFock> RefPFock;
typedef boost::shared_ptr<PTildeX> RefPTildeX;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<PMatrix1e> RefPMatrix1e;

using namespace std;
using namespace boost;

PSCF::PSCF(const RefPGeometry g) : geom_(g), overlap_(new POverlap(g)), hcore_(new PHcore(g)) {

  RefPTildeX tildex_tmp(new PTildeX(overlap_)); 
  tildex_ = tildex_tmp;

  direct_ = false;
  eig_ = new double[(2 * geom_->K() + 1) * geom_->nbasis()];

}

PSCF::~PSCF() {
  delete[] eig_;
}

void PSCF::compute() {
  init_schwarz();

  const string indent = "  ";
  const string space3 = "   ";

  RefPFock previous_fock;
  {
    RefPFock fock(new PFock(geom_, hcore_)); 
    previous_fock = fock;

    PMatrix1e intermediate = *tildex_ % *fock * *tildex_;
    intermediate.hermite();
    intermediate.diagonalize(eig_);
    RefPCoeff new_coeff(new PCoeff(*tildex_ * intermediate));
    RefPMatrix1e den(new PMatrix1e(new_coeff->form_density_rhf()));
    aodensity_ = den;
  }

  RefPMatrix1e densitychange = aodensity_; // assumes hcore guess...

  PDIIS<PMatrix1e> diis(5, 0.0);
  cout << indent << "=== Periodic RHF iteration (" + geom_->basisfile() + ")===" << endl << indent << endl;

  for (int iter = 0; iter != MAX_ITER_SCF; ++iter) { 
    int start = ::clock();

    RefPFock fock(new PFock(geom_, previous_fock, densitychange, schwarz_, geom_->S(), direct_, ao_eri_));
    previous_fock = fock;

    PMatrix1e intermediate = *tildex_ % *fock * *tildex_;
    intermediate.hermite();
    intermediate.diagonalize(eig_);
    RefPCoeff new_coeff(new PCoeff(*tildex_ * intermediate));
    coeff_ = new_coeff;

    RefPMatrix1e new_density(new PMatrix1e(coeff_->form_density_rhf()));
    new_density->real();
    RefPMatrix1e error_vector;
/*  if (iter < 10) { 
      PMatrix1e denft = aodensity_->ft();
      PMatrix1e ovrft = overlap_->ft();
      RefPMatrix1e evr(new PMatrix1e(*fock * denft * ovrft - ovrft * denft * *fock));
      error_vector = evr;
    } else {
*/  {
      RefPMatrix1e evr(new PMatrix1e(*new_density - *aodensity_));
      error_vector = evr;
    }


    int end = ::clock();
    const double energy = obtain_energy(*hcore_, fock->bft(), *new_density); 
    const double error = error_vector->rms();

    cout << indent << setw(5) << iter << setw(18) << fixed << setprecision(10) << energy << space3
                                      << setw(15) << error << setw(15) << setprecision(2) << (end - start) * 1.0e-6 << endl; 

    if (error < SCF_THRESH) {
      cout << indent << endl << indent << "  * SCF iteration converged." << endl << endl;
      break;
    }

    RefPMatrix1e diis_density;
#define SKIP
#ifndef SKIP
    if (iter > 5 && iter % 2 == 1) {
#endif
      diis_density = diis.extrapolate(make_pair(new_density, error_vector));
#ifndef SKIP
    } else {
      const double a = 1.0;
      new_density->scale(a);
      RefPMatrix1e tmp(new PMatrix1e(*new_density + *aodensity_));
      tmp->scale(1.0 / (1.0 + a));
      diis_density = tmp;
    }
#endif
    RefPMatrix1e dtmp(new PMatrix1e(*diis_density - *aodensity_));
    densitychange = dtmp; 
    aodensity_ = diis_density;
  }

  print_eig(eig_);

}


void PSCF::init_schwarz() {
  const vector<RefAtom> atoms = geom_->atoms(); 
  vector<RefShell> basis; 
  for (vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter) {
    const vector<RefShell> tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
  }
  const int size = basis.size(); // the number of shells per unit cell

  schwarz_.resize(size * size * (2 * geom_->K() + 1));
  int mcount = 0;
  for (int m = - geom_->K(); m <= geom_->K(); ++m, ++mcount) { 
    const double disp[3] = {0.0, 0.0, m * geom_->A()};
    for (int i0 = 0; i0 != size; ++i0) { // center unit cell
      const RefShell b0 = basis[i0];
      for (int i1 = 0; i1 != size; ++i1) {
        const RefShell b1 = basis[i1]->move_atom(disp);

        vector<RefShell> input;
        input.push_back(b0);
        input.push_back(b1);
        input.push_back(b0);
        input.push_back(b1);
        ERIBatch eribatch(input, 0.0);
        eribatch.compute();
        const double* eridata = eribatch.data();
        const int datasize = eribatch.data_size();
        double cmax = 0.0;
        for (int xi = 0; xi != datasize; ++xi, ++eridata) {
          const double absed = fabs(*eridata);
          if (absed > cmax) cmax = absed;
        }
        schwarz_[mcount * size * size + i0 * size + i1] = cmax;
      }
    }
  }
}


void PSCF::print_eig(const double* eig) const {
  for (int i = 0; i >= -geom_->K(); i -= 10) {
    const int ii = i + geom_->K();
    const int outsize = min(10, geom_->K()+1+i);
    const int bsize = min(geom_->nbasis(), 15);
    if (i != 0) cout << endl;
    cout << "      K ";
    for (int j = 0; j != outsize; ++j) {
      cout << setw(12) << abs(i)+j;
    }
    cout << endl;
    for (int k = 0; k != bsize; ++k) {
      if (k == geom_->nocc()/2 - 1) {
        cout << "   HOCO ";
      } else if (k == geom_->nocc()/2) {
        cout << "   LUCO ";
      } else {
        cout << "        ";
      }
      for (int j = 0; j != outsize; ++j) {
        cout << fixed << setw(12) << setprecision(6) << eig[(ii-j) * geom_->nbasis() + k];
      }
      cout << endl;
    }
  }
  cout << endl;
}


const double PSCF::obtain_energy(const PMatrix1e& hcore, const PMatrix1e& fock, const PMatrix1e& density) {

  const complex<double>* dat0 = hcore.data()->front();
  const complex<double>* dat1 = fock.data()->front();
  const complex<double>* dat2 = density.data()->front();
  
  double en = 0.0;
  const int s = geom_->S();
  for (int i = -s; i <= s; ++i) {  
    assert(hcore.blocksize() == fock.blocksize() && fock.blocksize() == density.blocksize());
    const complex<double>* cdat0 = dat0 + (i + s) * hcore.blocksize();
    const complex<double>* cdat1 = dat1 + (i + s) * fock.blocksize();
    const complex<double>* cdat2 = dat2 + (-i + s) * density.blocksize();
    const int n = fock.nbasis();
    for (int j = 0, cnt1 = 0; j != n; ++j) {
      for (int k = 0, cnt2 = j; k != n; ++k, ++cnt1, cnt2 += n) {
        en += cdat2[cnt2].real() * (cdat0[cnt1] + cdat1[cnt1]).real();
      }
    }
  }

  en += geom_->nuclear_repulsion();

  return en;

}

