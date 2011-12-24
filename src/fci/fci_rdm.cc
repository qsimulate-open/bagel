//
// Author: Toru Shiozaki
// Date  : Dec 2011
//

#include <src/fci/fci.h>
#include <src/fci/rdm.h>

using namespace std;
static const int unit = 1;
static const double one = 1.0;
static const double zero = 0.0;

void FCI::compute_rdm12(const int ist) {
  shared_ptr<Civec> cc = cc_->data(ist);
  
  // we need expanded lists
  const_phis_<0>(stringa_, phia_, false);
  const_phis_<1>(stringb_, phib_, false);

  const int la = cc->lena();
  const int lb = cc->lenb();
  const int len = la * lb;
  const int ij = norb_ * norb_;

  // creating new scratch dir.
  shared_ptr<Dvec> d(new Dvec(lb, la, ij));
  d->zero();
  sigma_2a1(cc, d);
  sigma_2a2(cc, d);

  // 1RDM
  shared_ptr<RDM1> rdm1(new RDM1(norb_));
  dgemv_("T", &len, &ij, &one, d->data(0)->first(), &len, cc->first(), &unit, &zero, rdm1->first(), &unit); 
  // 2RDM
  shared_ptr<RDM2> rdm2(new RDM2(norb_));
  dgemm_("T", "N", &ij, &ij, &len, &one, d->data(0)->first(), &len, d->data(0)->first(), &len,
          &zero, rdm2->first(), &ij);
  // put int diagonal into 2RDM
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      for (int j = 0; j != norb_; ++j) {
        rdm2->element(j,k,k,i) -= rdm1->element(j,i); 
      }
    }
  }

#if 0
  rdm1->print();
  rdm2->print();
#endif
#if 1
  // recomputing energy
  const int mm = norb_*norb_;
  const int nn = mm*mm;
  cout << endl << "     recomputing energy using RDMs : " << setprecision(12) << setw(18) << 
            geom_->nuclear_repulsion() + ddot_(&mm, jop_->mo1e_unpacked_ptr(), &unit, rdm1->first(), &unit)
                                       + 0.5*ddot_(&nn, jop_->mo2e_unpacked_ptr(), &unit, rdm2->first(), &unit) << endl; 
#endif
}

