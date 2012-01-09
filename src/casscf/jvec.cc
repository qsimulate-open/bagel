//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <src/casscf/jvec.h>

using namespace std;

Jvec::Jvec(shared_ptr<FCI> fci, shared_ptr<Coeff> coeff, const size_t nclosed, const size_t nact, const size_t nvirt) {

  shared_ptr<DensityFit> df = fci->geom()->df();
  const double* const cc = coeff->data();
  const size_t nocc = nclosed+nact;

  shared_ptr<RDM<1> > rdm1_av = fci->rdm1_av();
  shared_ptr<RDM<2> > rdm = fci->rdm2_av();

  // half_ = J^{-1/2}(D|ix)
  half_ = df->compute_half_transform(cc, nocc)->apply_J();

  // jvec_ = J^{-1/2} (D|ij) Gamma_ij,kl
  {
    shared_ptr<DF_Full> in = half_->compute_second_transform(cc,nocc);

    // for the time being, I form the entire 2RDM 
    unique_ptr<double[]> rdm2all(new double[nocc*nocc*nocc*nocc]); 
    fill(rdm2all.get(), rdm2all.get()+nocc*nocc*nocc*nocc, 0.0);
    {
      // closed-closed
      for (int i = 0; i != nclosed; ++i) {
        for (int j = 0; j != nclosed; ++j) {
          rdm2all[j+nocc*(j+nocc*(i+nocc*i))] += 4.0;
          rdm2all[j+nocc*(i+nocc*(i+nocc*j))] -= 2.0;
        }
      }
      // active-active
      for (int i = 0; i != nact; ++i) {
        for (int j = 0; j != nact; ++j) {
          for (int k = 0; k != nact; ++k) {
            dcopy_(nact, rdm->data()+nact*(k+nact*(j+nact*i)),1, rdm2all.get()+nclosed+nocc*(k+nclosed+nocc*(j+nclosed+nocc*(i+nclosed))),1);
          }
        }
      }
      // active-closed
      for (int i = 0; i != nclosed; ++i) {
        for (int j = 0; j != nact; ++j) {
          for (int k = 0; k != nact; ++k) {
            rdm2all[i+nocc*(i+nocc*(k+nclosed+nocc*(j+nclosed)))] += rdm1_av->element(k,j) * 2.0;
            rdm2all[k+nclosed+nocc*(j+nclosed+nocc*(i+nocc*(i)))] += rdm1_av->element(k,j) * 2.0;
            rdm2all[i+nocc*(k+nclosed+nocc*(j+nclosed+nocc*(i)))] -= rdm1_av->element(k,j);
            rdm2all[k+nclosed+nocc*(i+nocc*(i+nocc*(j+nclosed)))] -= rdm1_av->element(k,j);
          }
        }
      }
    }
    jvec_ = in->apply_2rdm(rdm2all.get())->apply_J();
    rdm2_all_ = move(rdm2all);
#if 0
double* d = rdm2_all_.get();
for (int i = 0; i != nocc; ++i) {
for (int j = 0; j != nocc; ++j) {
for (int k = 0; k != nocc; ++k) {
for (int l = 0; l != nocc; ++l, ++d) {
if (abs(*d) > 1.0e-10) cout << i << " " << j << " " << k << " " << l << " " << setprecision(10) << *d << endl;
}}}}
#endif
  }

};
