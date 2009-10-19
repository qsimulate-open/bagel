//
// Author : Toru Shiozaki
// Date   : August 2009
//

#include <cstring>
#include <iostream>
#include <algorithm>
#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/macros.h>
#include <src/scf/scf_macros.h>
#include <src/slater/slaterbatch.h>
#include <src/util/paircompfile.h>
#include <src/util/pcompcabsfile.h>
#include <src/util/pmofile.h>

typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<PGeometry> RefGeom;
typedef boost::shared_ptr<PCoeff> RefPCoeff;
typedef boost::shared_ptr<Shell> RefShell;

using namespace std;
using namespace boost;

PMP2::PMP2(const RefGeom g, const RefPCoeff co, const double* eg, const shared_ptr<PCompFile<ERIBatch> > fl)
 : geom_(g), coeff_(co), eig_(eg), ao_eri_(fl) {

  nfrc_ = geom_->nfrc() / 2;
  nocc_ = geom_->nocc() / 2;
  nocc_act_ = nocc_ - nfrc_;
  nbasis_ = geom_->nbasis();
  nvir_ = nbasis_ - nocc_;
  noovv_ = nocc_act_ * nocc_act_ * nbasis_ * nbasis_;
  ncabs_ = geom_->ncabs();

  assert(geom_->ncabs() != 0);
// just to check
//const vector<RefAtom> tmp = geom_->cabs_atoms();
//for (vector<RefAtom>::const_iterator iter = tmp.begin(); iter != tmp.end(); ++iter) (*iter)->print_basis();

}


PMP2::~PMP2() {

}


void PMP2::compute() {

  // AO ERI has been computed in the SCF class.

  cout << "  === Periodic MP2 calculation ===" << endl << endl;
  // Fully transform aa/ii integrals and dump them to disk (... forcus is on MP2-R12).
  eri_ii_pp_ = ao_eri_->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                             0, nbasis_, 0, nbasis_, "ERI (pp/ii)");
  eri_ii_pp_->sort_inside_blocks();

  // Compute the conventional MP2 contribution
  compute_conv_mp2();

  cout << "  === Periodic MP2-R12 calculation ===" << endl << endl;
  // Calculate Yukawa potential integrals
  shared_ptr<PairCompFile<SlaterBatch> > stg_yp(new PairCompFile<SlaterBatch>(geom_, 1.5, "Slater and Yukawa ints"));
  stg_yp->store_integrals();
  stg_yp->reopen_with_inout();
  shared_ptr<PCompFile<SlaterBatch> > stg = stg_yp->first();
  shared_ptr<PCompFile<SlaterBatch> > yp  = stg_yp->second();

  // V intermediate OBS part
  typedef shared_ptr<PMOFile<complex<double> > > RefPMOFile;
  RefPMOFile yp_ii_ii = yp->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                                 nfrc_, nocc_, nfrc_, nocc_, "Yukawa (ii/ii)");
  yp_ii_ii->sort_inside_blocks();
  RefPMOFile stg_ii_pp = stg->mo_transform(coeff_, nfrc_, nocc_, nfrc_, nocc_,
                                                   0, nbasis_, 0, nbasis_, "Slater (pp/ii)");
  stg_ii_pp->sort_inside_blocks();

  RefPMOFile stg_dag_times_eri_ii_ii = stg_ii_pp->contract(eri_ii_pp_, nfrc_, nocc_,
                                                                       nfrc_, nocc_, "F * v (ii/ii)");
  RefPMOFile V_obs(new PMOFile<complex<double> >(*yp_ii_ii - *stg_dag_times_eri_ii_ii)); 

//V_obs->print();
//yp_ii_ii->print();

  complex<double> en_vt = V_obs->get_energy_one_amp();
  // Direct contribution to R12 energy
  cout << "  F12 energy (Vt): " << setprecision(10) << en_vt.real() << endl << endl;

  // CABS integrals
  shared_ptr<PCompCABSFile<ERIBatch> >eri_cabs_(new PCompCABSFile<ERIBatch>(geom_, false, "ERI CABS"));
  eri_cabs_->store_integrals();
  eri_cabs_->reopen_with_inout();

  // Construction of CABS;
  shared_ptr<PMatrix1e> cabs_obs;
  shared_ptr<PMatrix1e> cabs_aux;
  {
    typedef shared_ptr<PMatrix1e> RefMatrix;

    // Form RI space which is a union of OBS and CABS.
    RefGeom union_geom(new PGeometry(*geom_));
    union_geom->merge_obs_cabs();

    shared_ptr<POverlap> union_overlap(new POverlap(union_geom));
    shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));
    RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim(), coeff_->mdim()));

    // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
    RefMatrix tmp(new PMatrix1e(*ri_coeff % union_overlap->ft() * *ri_reshaped));
    RefMatrix U(new PMatrix1e(geom_, tmp->ndim(), tmp->ndim()));
    RefMatrix V(new PMatrix1e(geom_, tmp->mdim(), tmp->mdim()));
    tmp->svd(U, V);

    RefMatrix Ured(new PMatrix1e(U, tmp->mdim()));
    RefMatrix cabs_coeff(new PMatrix1e(*ri_coeff * *Ured));

    pair<RefMatrix, RefMatrix> cabs_coeff_spl = cabs_coeff->split(geom_->nbasis(), geom_->ncabs());
    cabs_coeff->print();
    cabs_obs = cabs_coeff_spl.first;
    cabs_aux = cabs_coeff_spl.second;

    // Setting actual number of CABS.
    ncabs_ = cabs_obs->mdim();
  }
  assert(cabs_obs->mdim() == ncabs_);
  RefPMOFile eri_ii_ii = ao_eri_->mo_transform_cabs_obs(coeff_, cabs_obs,
                                                        nfrc_, nocc_, nfrc_, nocc_,
                                                        nfrc_, nocc_, 0, ncabs_, "V^ia'_ii, OBS part");


}


void PMP2::compute_conv_mp2() {
  cout << "  Now computing mp2 energy" << endl << endl;
  const int K = eri_ii_pp_->K();
  const int maxK1 = max(K, 1);
  const int KK = K + K;

  complex<double>* workoovv1 = new complex<double>[noovv_];
  complex<double>* workoovv2 = new complex<double>[noovv_];

  double energy = 0.0;
  int nbja = 0;
  for (int kb = -K; kb < maxK1; ++kb) {
    const double* eb = eig_ + (kb + K) * nbasis_;

    for (int kj = -K; kj != maxK1; ++kj) {
      const double* ej = eig_ + (kj + K) * nbasis_;

      for (int ka = -K; ka != maxK1; ++ka, ++nbja) {
        const double* ea = eig_ + (ka + K) * nbasis_;

        int ki = ka + kb - kj; 
        if (ki < -K) ki += K * 2; 
        else if (ki >= K) ki -= K * 2; 
        const double* ei = eig_ + (ki + K) * nbasis_;

        assert(K == 0 || (ki < K && ki >= -K));

        {
          const long tmp2 = kb + K + KK * (kj + K + KK * (ka + K)); 
          eri_ii_pp_->get_block(noovv_ * nbja, noovv_, workoovv1);
          eri_ii_pp_->get_block(noovv_ * tmp2, noovv_, workoovv2);
        }

        #pragma omp parallel for reduction(+:energy) 
        for (int i = nfrc_; i < nocc_; ++i) {
          const int xi = i - nfrc_;
          for (int j = nfrc_, xj = 0; j != nocc_; ++j, ++xj) {
            const int ijoffset = nbasis_ * (xj + nocc_act_ * xi);
            for (int a = nocc_; a != nbasis_; ++a) {
              for (int b = nocc_; b != nbasis_; ++b) {
                int cnt1 = b + nbasis_ * (a + ijoffset);
                int cnt2 = a + nbasis_ * (b + ijoffset);
                const complex<double> v1 = workoovv1[cnt1];
                const complex<double> v2 = workoovv2[cnt2];
                energy += (real(v1 * conj(v1 + v1 - v2))) / (ei[i] + ej[j] - ea[a] - eb[b]);
              }
            }
          }
        }
      }
    }
  }
  cout << "  MP2 energy: " << setprecision(10) << energy / (max(pow(2.0 * K, 3.0), 1.0))  << endl << endl;
  delete[] workoovv1;
  delete[] workoovv2;

}

