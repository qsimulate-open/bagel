
#include <src/asd/orbopt/ras/bfgs.h>

using namespace std;
using namespace bagel;

// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
// G_(i<a)
void ASDRASBFGS::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<ASDRASRotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}

// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
// G_(t<a)
void ASDRASBFGS::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASDRASRotFile> sigma) const {
  if (!nvirt_ || !nact_) return;
  dgemm_("N", "T", nvirt_, nact_, nact_, 2.0, cfock->element_ptr(nocc_,nclosed_), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_);
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
//  daxpy_(nvirt_, 2.0*occup_[i], cfock->element_ptr(nocc_, i+nclosed_), 1, target, 1);
    daxpy_(nvirt_, 2.0, qxr->element_ptr(nocc_, i), 1, target, 1);
  }
}

// grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
// G_(i<t)
void ASDRASBFGS::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASDRASRotFile> sigma) const {
  if (!nclosed_ || !nact_) return;
  {
    double* target = sigma->ptr_ca();
    for (int i = 0; i != nact_; ++i, target += nclosed_) {
    //daxpy_(nclosed_, 4.0-2.0*occup_[i], cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, afock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, -2.0, qxr->element_ptr(0, i), 1, target, 1);
    }
    //-2 cfock_iu * D_ur
    dgemm_("T", "N", nclosed_, nact_, nact_, -2.0, cfock->element_ptr(nclosed_,0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, sigma->ptr_ca(), nclosed_);
  }
}

// grad(t/t)
// G_(t(A)<t(B))
void ASDRASBFGS::grad_aa(shared_ptr<const Matrix> mcfock, shared_ptr<ASDRASRotFile> sigma) const {
  if (!nact_) return;
  double* target = sigma->ptr_aa();
  for (int jb = nactA_; jb != nact_; ++jb) { //B
    for (int ia = 0; ia != nactA_; ++ia, ++target) { //A
      *target = -2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
      //cout << "Grad(" << i << "," << j << ") = " << 2.0*(mcfock->element(ia,jb) - mcfock->element(ia,jb)) << endl;
    }
  }
}

// grad(t/t) /with previous gradient as constant
void ASDRASBFGS::grad_aa_with_preg(shared_ptr<const Matrix> mcfock, shared_ptr<ASDRASRotFile> sigma) const {
  if (!nact_) return;
  double* target = sigma->ptr_aa();
  double* preg_p = preg_->data();
  for (int i = 0; i != nactB_; ++i) { //B
    for (int j = 0; j != nactA_; ++j, ++target, ++preg_p) { //A
      *target = *preg_p - 2.0*(mcfock->element(j,i+nactA_) - mcfock->element(i+nactA_,j));
      *preg_p = 2.0*(mcfock->element(j,i+nactA_) - mcfock->element(i+nactA_,j));
    //*target = 2.0*(mcfock->element(j,i+nactA_) - mcfock->element(i+nactA_,j));
      cout << "Grad(" << i << "," << j << ") = " << 2.0*(mcfock->element(j,i+nactA_) - mcfock->element(i+nactA_,j)) << endl;
    //*target = - mcfock->element(j,i+nactA_) + mcfock->element(i+nactA_,j);
    }
  }
}
