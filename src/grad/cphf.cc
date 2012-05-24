//
// Author : Toru Shiozaki
// Date   : May 2012
//

#include <src/grad/cphf.h>
#include <cassert>

#define CPHF_MAX_ITER 100
#define CPHF_THRESH 1.0e-8

using namespace std;

CPHF::CPHF(const shared_ptr<const Matrix1e> grad, const vector<double>& eig, const shared_ptr<DF_Half> h,
           const shared_ptr<const Reference> r, const int n)
: solver_(new Linear<Matrix1e>(CPHF_MAX_ITER, grad)), grad_(grad), eig_(eig), halfjj_(h), ref_(r), geom_(r->geom()), ncore_(n) {

}


shared_ptr<Matrix1e> CPHF::solve() const {

  const size_t naux = geom_->naux();
  const size_t nocc = geom_->nele() / 2 - ncore_;
  if (nocc < 1) throw runtime_error("no correlated electrons"); 
  const size_t nvirt = geom_->nbasis() - nocc - ncore_;
  if (nvirt < 1) throw runtime_error("no virtuals orbitals"); 
  assert(geom_->nbasis() == ref_->coeff()->mdim());

  const size_t nbasis = geom_->nbasis();

  const double* const coeff = ref_->coeff()->data() + ncore_*nbasis;
  const double* const vcoeff = coeff + nocc*nbasis;

  shared_ptr<Matrix1e> t(new Matrix1e(geom_));
  for (int i = 0; i != nocc; ++i)
    for (int a = nocc; a != nvirt+nocc; ++a)
      t->element(a,i) = grad_->element(a,i) / (eig_[a]-eig_[i]);

  unique_ptr<double[]> jri(new double[nbasis*nocc]);
  unique_ptr<double[]> jai(new double[nvirt*nocc]);
  unique_ptr<double[]> kia(new double[nvirt*nocc]);

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
    solver_->orthog(t);

    shared_ptr<Matrix1e> sigma(new Matrix1e(geom_));
    // one electron part
    sigma->zero();
    for (int i = 0; i != nocc; ++i)
      for (int a = nocc; a != nvirt+nocc; ++a)
        sigma->element(a,i) = (eig_[a]-eig_[i]) * t->element(a,i);

    // J part
    shared_ptr<Matrix1e> pbmao(new Matrix1e(geom_));
    unique_ptr<double[]> pms(new double[nbasis*nocc]);
    dgemm_("T", "T", nocc, nbasis, nvirt, 1.0, t->data()+nocc, nbasis, vcoeff, nbasis, 0.0, pms.get(), nocc);
    dgemm_("N", "N", nbasis, nbasis, nocc, 1.0, coeff, nbasis, pms.get(), nocc, 0.0, pbmao->data(), nbasis);
    pbmao->symmetrize();

    unique_ptr<double[]> jrs = geom_->df()->compute_Jop(pbmao->data());
    dgemm_("N", "N", nbasis, nocc, nbasis, 1.0, jrs.get(), nbasis, coeff, nbasis, 0.0, jri.get(), nbasis);
    dgemm_("T", "N", nvirt, nocc, nbasis, 4.0, vcoeff, nbasis, jri.get(), nbasis, 0.0, jai.get(), nvirt); 
    // K part
    unique_ptr<double[]> kir = halfjj_->compute_Kop_1occ(pbmao->data());
    dgemm_("N", "N", nocc, nvirt, nbasis, -2.0, kir.get(), nocc, vcoeff, nbasis, 0.0, kia.get(), nocc); 
    for (int i = 0; i != nocc; ++i)
      for (int a = 0; a != nvirt; ++a)
        sigma->element(a+nocc,i) += jai[a+nvirt*i] + kia[i+nocc*a];

    t = solver_->compute_residual(t, sigma);

    // TODO to be controlled by the input
    if (t->norm() < CPHF_THRESH) break;

    for (int i = 0; i != nocc; ++i)
      for (int a = nocc; a != nvirt+nocc; ++a)
        t->element(a,i) /= (eig_[a]-eig_[i]);

  }
  t = solver_->civec();
  t->fill_upper();
  return t;

}
