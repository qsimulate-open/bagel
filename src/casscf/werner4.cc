//
// Author : Toru Shiozaki
// Date   : April 2012
//

// implements second-order Newton-Raphson CASSCF

#include <src/casscf/werner4.h>

using namespace std;


void Werner4::compute() {
  const string indent = "  ";
  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

  // initializing Hcore matrix
  shared_ptr<Hcore> hc(new Hcore(geom_));
  shared_ptr<Fock<1> > fc(new Fock<1>(geom_, hc));
  hcore_ = fc;

  // macro iteration starts here
  for (int iter = 0; iter != max_iter_; ++iter) {
    int start = ::clock();
    // performs FCI, which also computes one-index transformed integrals
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    vector<double> energy = fci_->energy();

    // start with U=1
    shared_ptr<Matrix1e> U(new Matrix1e(geom_));
    U->unit();

    // compute jvec
    shared_ptr<Jvec> jvec(new Jvec(fci_, ref_->coeff(), nclosed_, nact_, nvirt_));

    // compute the G operator
    // TODO for the time being, it is done in a very robust way // note the redundancy (with Jvec above)
    JKop jk(geom_->df(), ref_->coeff(), hcore_, fci_, nocc_, nclosed_, nact_);

    // micro iteration starts here
    int miter;
    double error = 0;
    for (miter = 0; miter != max_micro_iter_; ++miter) {

      // compute the (initial) B vector.
      shared_ptr<Matrix1e> bvec = compute_bvec(jvec, U, ref_->coeff());
      shared_ptr<Matrix1e> bvec0 = bvec;
      // update C = 1/2(A+A^dagger) = 1/2(U^dagger B + B^dagger U)
      shared_ptr<Matrix1e> C(new Matrix1e((*U % *bvec + *bvec % *U)*0.5));

      // for debug
      shared_ptr<Matrix1e> t(new Matrix1e(*U));
      for (int i = 0; i != t->ndim(); ++i) t->element(i,i) -= 1.0;
      const double energy_micro = (bvec0->ddot(*t)*0.5+ bvec->ddot(*t)*0.5)+ energy[0];
      cout << "micro energy " << setw(20) << setprecision(10) << energy_micro << " " << bvec0->ddot(*t) << " " << bvec->ddot(*t) << endl;

      // compute gradient
      shared_ptr<Matrix1e> grad(new Matrix1e(*U%*bvec-*bvec%*U));
      grad->purify_redrotation(nclosed_,nact_,nvirt_);
      const double error_micro = grad->ddot(*grad)/grad->size();
      if (miter == 0) error = error_micro;
      cout << "   -- iter " << setw(4) << miter << "  residual norm : " << setw(20) << setprecision(10) << error_micro << endl;

      // solving Eq. 30 of Werner and Knowles
      // fixed U, solve for dR.
      for (int mmiter = 0; mmiter != max_mmicro_iter_; ++mmiter) {

      }
break;
    }

break;
  }
};
