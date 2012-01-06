//
// Author : Toru Shiozaki
// Date   : Jan 2012
//


#include <src/casscf/werner.h>
#include <src/casscf/jvec.h>
#include <src/scf/matrix1e.h>
#include <src/casscf/rotfile.h>
#include <src/util/davidson.h>

using namespace std;

#define DF 1

void WernerKnowles::compute() {
  const string indent = "  ";

  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

  // initializing Hcore matrix (redundant copy, but I can live with it).
  shared_ptr<Matrix1e> hcore_;
  {
    shared_ptr<Hcore> hc(new Hcore(geom_));
    shared_ptr<Fock<DF> > fc(new Fock<DF>(geom_, hc));
    shared_ptr<Matrix1e> hcore(new Matrix1e(*ref_->coeff() % *fc * *ref_->coeff()));
    hcore_ = hcore;
  }

  // All the equation numbers refer to as those in Werner and Knowles, J. Chem. Phys. 82, 5053 (1985).
  // macro iteration
  for (int iter = 0; iter != 1; ++iter) {
//for (int iter = 0; iter != max_iter_; ++iter) {

    // performs FCI, which also computes one-index transformed integrals 
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    vector<double> energy = fci_->energy();

    shared_ptr<Jvec> jvec(new Jvec(fci_, ref_->coeff(), nclosed_, nact_, nvirt_));
    shared_ptr<RotFile> R(new RotFile(nclosed_, nact_, nvirt_));

    // initially U is unit
    shared_ptr<Matrix1e> U(new Matrix1e(geom_));
    U->unit();

    // unit matrix is also prepared for convenience
    shared_ptr<Matrix1e> one(new Matrix1e(*U));

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      shared_ptr<Matrix1e> T(new Matrix1e(*U - *one));
      // compute initial B (Eq. 19)
      shared_ptr<Matrix1e> bvec = compute_bvec(hcore_, fci_, jvec, U, T, ref_->coeff());

      // update C = 1/2(A+A^dagger) = 1/2(U^dagger B + B^dagger U)

      // initializing a Davidson manager
      DavidsonDiag<RotFile> davidson(1,max_micro_iter_); 

break;

      // solving Eq. 30 of Werner and Knowles
      for (int mmiter = 0; mmiter != max_micro_iter_; +mmiter) {
        if (mmiter != 0) {
          // update B

          // compute  C dR + dR C

        }

        // compute U^dagger B - B^dagger U

        // computes Epsilon

        // computes Eq.30

      }

    }

  }

}


// compute B according to Eq. (19).
// B = 2[h_rs U_sj D_ji + (rs|D)Usj <D|ji> + 2(rk|D)(D|ls)T_sj D_jl,ik]
shared_ptr<Matrix1e> WernerKnowles::compute_bvec(shared_ptr<Matrix1e> hcore, shared_ptr<FCI> fci, shared_ptr<Jvec> jvec,
                                                 shared_ptr<Matrix1e> u, shared_ptr<Matrix1e> t, shared_ptr<Coeff> cc) {
  shared_ptr<DensityFit> df = geom_->df();
  const int naux = df->naux();
  const int nbas = df->nbasis();
  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  {
    // first term
    Matrix1e all1(geom_);
    for (int i = 0; i != nclosed_; ++i) all1.element(i,i) = 2.0;
    for (int i = 0; i != nact_; ++i) dcopy_(nact_, fci->rdm1_av()->data()+nact_*i, 1, all1.element_ptr(nclosed_, i+nclosed_), 1); 
    shared_ptr<Matrix1e> buf(new Matrix1e(geom_));
    dgemm_("N", "N", nbasis_, nocc_, nbasis_, 1.0, hcore->data(), nbasis_, u->data(), nbasis_, 0.0, buf->data(), nbasis_);  
    dgemm_("N", "N", nbasis_, nocc_, nocc_, 2.0, buf->data(), nbasis_, all1.data(), nbasis_, 0.0, out->data(), nbasis_);
  }

  unique_ptr<double[]> tmp(new double[nbas*nbas]);
  // second term
  {
    Matrix1e Umat(*cc * *u);
    shared_ptr<DF_Half> half = df->compute_half_transform(Umat.data(), nocc_);
    half->form_2index(tmp, jvec->jvec()->apply_J(), 2.0, 0.0);
  }

  // thrid term
  if (t->norm() > 1.0e-15) {
    Matrix1e Tmat(*cc * *t);
    shared_ptr<DF_Full> full = jvec->half()->compute_second_transform(Tmat.data(), nocc_)->apply_2rdm(jvec->rdm2_all());
    jvec->half()->form_2index(tmp, full, 4.0, 1.0);
  }

  dgemm_("T", "N", nbasis_, nocc_, nbas, 1.0, cc->data(), nbas, tmp.get(), nbas, 1.0, out->data(), nbasis_);
  return out;
}

