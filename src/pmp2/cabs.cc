/*
 * cabs.cc
 *
 *  Created on: Oct 22, 2009
 *      Author: shiozaki
 */

#include <src/pmp2/pmp2.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/ptildex.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pcoeff.h>
#include <src/util/pcompcabsfile.h>

using namespace std;
using namespace boost;

//#define LOCAL_DEBUG

typedef shared_ptr<PMatrix1e> RefMatrix;
typedef shared_ptr<PGeometry> RefGeom;
typedef shared_ptr<PHcore> RefHcore;
typedef shared_ptr<PCoeff> RefCoeff;
typedef shared_ptr<PMOFile<complex<double> > > RefMOFile;

pair<RefCoeff, RefCoeff> PMP2::generate_CABS() {

  // Form RI space which is a union of OBS and CABS.
  RefGeom newgeom(new PGeometry(*geom_));
  union_geom_ = newgeom;
  union_geom_->merge_obs_cabs();

  shared_ptr<POverlap> union_overlap(new POverlap(union_geom_));
  shared_ptr<PTildeX> ri_coeff(new PTildeX(union_overlap));
  RefMatrix ri_reshaped(new PMatrix1e(coeff_, ri_coeff->ndim()));

  // SVD to project out OBS component. Note singular values are all 1 as OBS is a subset of RI space.
  RefMatrix tmp(new PMatrix1e(*ri_coeff % union_overlap->ft() * *ri_reshaped));

  const int tmdim = tmp->mdim();
  const int tndim = tmp->ndim();

  RefMatrix U(new PMatrix1e(geom_, tndim, tndim));
  RefMatrix V(new PMatrix1e(geom_, tmdim, tmdim));
  tmp->svd(U, V);

  RefMatrix Ured(new PMatrix1e(U, make_pair(tmdim, tndim)));
  RefCoeff coeff_cabs(new PCoeff(*ri_coeff * *Ured));
  coeff_cabs_ = coeff_cabs;

  RefMatrix coeff_fit(new PMatrix1e(coeff_, tndim));
  RefMatrix coeff_entire(new PMatrix1e(coeff_fit, coeff_cabs_));
  coeff_entire_ = coeff_entire;

  pair<RefCoeff, RefCoeff> cabs_coeff_spl = coeff_cabs_->split(geom_->nbasis(), geom_->ncabs());

  return cabs_coeff_spl;
}


const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_hJ() const {

#define USE_BUILDER
  // hcore contribution.
  RefHcore uhc(new PHcore(union_geom_, true));

#ifdef USE_BUILDER

  RefMatrix ao_hJ(new PMatrix1e((*uhc + *coulomb_runtime()).ft()));
  RefMatrix hJ(new PMatrix1e(*coeff_entire_ % *ao_hJ * *coeff_entire_));

#else
  RefMatrix aohcore(new PMatrix1e(uhc->ft()));
  RefMatrix mohcore(new PMatrix1e(*coeff_entire_ % *aohcore * *coeff_entire_));
  // coulomb obs_obs block
  RefMatrix coulomb_oo;
  {
    RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                                 0, nbasis_, 0, nocc_,
                                                 0, nbasis_, 0, nocc_, "h+J builder (OBS-OBS; pp)");


    coulomb_oo = eri_pI_pI->contract_density_J();
  }

  // coulomb obs_cabs block
  RefMatrix coulomb_oc;
  {
    RefMOFile eri_pI_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                 0, nbasis_, 0, nocc_,
                                                 0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; pp)");
    RefMOFile eri_pI_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                           0, nbasis_, 0, nocc_,
                                                           0, ncabs_, 0, nocc_, "h+J builder (OBS-CABS; px)");
    RefMOFile eri_pI_AI(new PMOFile<complex<double> >(*eri_pI_xI + *eri_pI_pI));

    coulomb_oc = eri_pI_AI->contract_density_J();
  }


  const double gamma = geom_->gamma();
  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, false, false, false, "ERI CABS(j)"));
  eri_cabs_d->store_integrals();
  eri_cabs_d->reopen_with_inout();

  // coulomb cabs_obs block
  RefMatrix coulomb_co;
  {
    RefMOFile eri_AI_pI = eri_obs_->mo_transform(cabs_obs_, coeff_, coeff_, coeff_,
                                                 0, ncabs_, 0, nocc_,
                                                 0, nbasis_, 0, nocc_, "hJ builder (CABS-OBS; pp)");
    *eri_AI_pI += *(eri_cabs_d->mo_transform_cabs_aux(cabs_aux_, coeff_, coeff_, coeff_,
                                                      0, ncabs_, 0, nocc_,
                                                      0, nbasis_, 0, nocc_, "hJ builder (CABS-OBS; xp)"));
    coulomb_co = eri_AI_pI->contract_density_J();
  }

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, true, false, false, "ERI CABS(ja)"));
  eri_cabs_t->store_integrals();
  eri_cabs_t->reopen_with_inout();

  // coulomb cabs_cabs block
  RefMatrix coulomb_cc;
  {
    RefMOFile eri_AI_AI = eri_cabs_t->mo_transform_cabs_aux(cabs_aux_, coeff_, cabs_aux_, coeff_,
                                                            0, ncabs_, 0, nocc_,
                                                            0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; xx)");
    *eri_AI_AI += *(eri_obs_->mo_transform(cabs_obs_, coeff_, cabs_obs_, coeff_,
                                           0, ncabs_, 0, nocc_,
                                           0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; pp)"));
    *eri_AI_AI += *(eri_cabs_d->mo_transform_cabs_aux(cabs_aux_, coeff_, cabs_obs_, coeff_,
                                                      0, ncabs_, 0, nocc_,
                                                      0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; xp)"));
    *eri_AI_AI += *(eri_cabs_->mo_transform_cabs_aux(cabs_obs_, coeff_, cabs_aux_, coeff_,
                                                     0, ncabs_, 0, nocc_,
                                                     0, ncabs_, 0, nocc_, "hJ builder (CABS-CABS; px)"));
    coulomb_cc = eri_AI_AI->contract_density_J();
  }

  RefMatrix coulomb(new PMatrix1e(coulomb_oo->merge(coulomb_oc), coulomb_co->merge(coulomb_cc)));
  coulomb->conj();
  RefMatrix hJ(new PMatrix1e(*coulomb+*mohcore));
#endif

  RefMatrix h_hJ_o(new PMatrix1e(hJ, make_pair(0, geom_->nbasis())));
  RefMatrix h_hJ_c(new PMatrix1e(hJ, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_hJ_o_pair = h_hJ_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_hJ_c_pair = h_hJ_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_hJ_o_pair.first, h_hJ_o_pair.second, h_hJ_c_pair.first, h_hJ_c_pair.second);
}




const boost::tuple<RefMatrix, RefMatrix, RefMatrix, RefMatrix> PMP2::generate_K() const {
  // TODO INEFFICIENT CODE!!! Hartree matrix needs to be constructed in AO basis.

  // obs_obs block
  RefMatrix exchange_oo;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, coeff_, coeff_,
                                                 0, nocc_, 0, nbasis_,
                                                 0, nbasis_, 0, nocc_, "K builder (OBS-OBS; pp)");
    exchange_oo = eri_Ip_pI->contract_density_K();
  }

  // obs_cabs block
  RefMatrix exchange_oc;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, coeff_, cabs_obs_, coeff_,
                                                 0, nocc_, 0, nbasis_,
                                                 0, ncabs_, 0, nocc_, "K builder (OBS-CABS; pp)");
    RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, coeff_, cabs_aux_, coeff_,
                                                           0, nocc_, 0, nbasis_,
                                                           0, ncabs_, 0, nocc_, "K builder (OBS-CABS; px)");
    RefMOFile eri_Ip_AI(new PMOFile<complex<double> >(*eri_Ip_xI + *eri_Ip_pI));
    exchange_oc = eri_Ip_AI->contract_density_K();
  }

  const double gamma = geom_->gamma();
  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, false, false, false, "ERI CABS(j)"));
  eri_cabs_d->store_integrals();
  eri_cabs_d->reopen_with_inout();

  //  cabs_obs block
  RefMatrix exchange_co;
  {
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, coeff_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, nbasis_, 0, nocc_, "K builder (CABS-OBS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, coeff_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, nbasis_, 0, nocc_, "K builder (CABS-OBS; xp)");
    RefMOFile eri_IA_pI(new PMOFile<complex<double> >(*eri_Ix_pI + *eri_Ip_pI));
    exchange_co = eri_IA_pI->contract_density_K();
  }

  shared_ptr<PCompCABSFile<ERIBatch> >
    eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, false, true, true, false, false, "ERI CABS(ja)"));
  eri_cabs_t->store_integrals();
  eri_cabs_t->reopen_with_inout();

  // cabs-cabs block
  RefMatrix exchange_cc;
  {
    RefMOFile eri_Ix_xI = eri_cabs_t->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_aux_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-CABS; xx)");
    RefMOFile eri_Ip_pI = eri_obs_->mo_transform(coeff_, cabs_obs_, cabs_obs_, coeff_,
                                                 0, nocc_, 0, ncabs_,
                                                 0, ncabs_, 0, nocc_, "K builder (CABS-CABS; pp)");
    RefMOFile eri_Ix_pI = eri_cabs_d->mo_transform_cabs_aux(coeff_, cabs_aux_, cabs_obs_, coeff_,
                                                            0, nocc_, 0, ncabs_,
                                                            0, ncabs_, 0, nocc_, "K builder (CABS-CABS; xp)");
    RefMOFile eri_Ip_xI = eri_cabs_->mo_transform_cabs_aux(coeff_, cabs_obs_, cabs_aux_, coeff_,
                                                           0, nocc_, 0, ncabs_,
                                                           0, ncabs_, 0, nocc_, "K builder (CABS-CABS; px)");
    RefMOFile eri_IA_AI(new PMOFile<complex<double> >(*eri_Ip_pI + *eri_Ix_pI + *eri_Ip_xI + *eri_Ix_xI));
    exchange_cc = eri_IA_AI->contract_density_K();
  }

  // first form entire matrix
  RefMatrix exchange(new PMatrix1e(exchange_oo->merge(exchange_oc), exchange_co->merge(exchange_cc)));
  // this is an important step!!!
  exchange->conj();

  RefMatrix h_exchange_o(new PMatrix1e(exchange, make_pair(0, geom_->nbasis())));
  RefMatrix h_exchange_c(new PMatrix1e(exchange, make_pair(geom_->nbasis(), geom_->nbasis()+geom_->ncabs())));

  pair<RefMatrix, RefMatrix> h_exchange_o_pair = h_exchange_o->split(geom_->nbasis(), geom_->ncabs());
  pair<RefMatrix, RefMatrix> h_exchange_c_pair = h_exchange_c->split(geom_->nbasis(), geom_->ncabs());

  return make_tuple(h_exchange_o_pair.first, h_exchange_o_pair.second, h_exchange_c_pair.first, h_exchange_c_pair.second);
}



RefMatrix PMP2::coulomb_runtime_OBS() const {

  RefMatrix density(new PMatrix1e(coeff_->form_density_rhf()));
  const complex<double>* den = density->data()->front();

  RefMatrix coulomb_real_space(new PMatrix1e(geom_));
  complex<double>* data = coulomb_real_space->data()->front();

  size_t allocsize = eri_obs_->max_num_int();
  double* diskdata = new double[allocsize];

  const int s = geom_->S();
  const int l = geom_->L();

  long file_position = 0l;
  size_t mcnt = 0lu;
  for (int m1 = -s; m1 <= s; ++m1) {
    for (int m2 = 0; m2 <= l; ++m2) { // use bra-ket symmetry!!!
      for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++mcnt) {

        const int k = geom_->K();
        const size_t b = coulomb_real_space->blocksize();
        assert(coulomb_real_space->blocksize() == density->blocksize());
        const int n = geom_->nbasis();

        const int m1________k____b = (m1      + k) * b;
        const int m1___m2___k____b = (m1 - m2 + k) * b;
        const int m2___m1___k____b = (m2 - m1 + k) * b;
        const int m2___m3___k____b = (m2 - m3 + k) * b;
        const int m3___m2___k____b = (m3 - m2 + k) * b;
        const int m3________k____b = (m3      + k) * b;
        const int _____m1___k____b = (   - m1 + k) * b;
        const int _____m3___k____b = (   - m3 + k) * b;

        {
          eri_obs_->get_block(file_position, eri_obs_->num_int_each(mcnt), diskdata);
          file_position += eri_obs_->num_int_each(mcnt);
          const double* cdata = diskdata;

          const int size = eri_obs_->basissize(); // number of shells
          for (int i0 = 0; i0 != size; ++i0) {
            const int b0offset = eri_obs_->offset(i0);
            const int b0size = eri_obs_->nbasis(i0);

            for (int i1 = 0; i1 != size; ++i1) {
              const int b1offset = eri_obs_->offset(i1);
              const int b1size = eri_obs_->nbasis(i1);

              for (int i2 = 0; i2 != size; ++i2) {
                const int b2offset = eri_obs_->offset(i2);
                const int b2size = eri_obs_->nbasis(i2);

                for (int i3 = 0; i3 != size; ++i3) {
                  const int b3offset = eri_obs_->offset(i3);
                  const int b3size = eri_obs_->nbasis(i3);

                  const double integral_bound = eri_obs_->schwarz(((m1      + k) * size + i0) * size + i1)
                                              * eri_obs_->schwarz(((m3 - m2 + k) * size + i2) * size + i3);
                  const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                  if (skip_schwarz) continue;

                  if (m2 != 0) {
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell
                      for (int j1 = b1offset, j1n = b1offset * n; j1 != b1offset + b1size; ++j1, j1n += n) {
                        for (int j2 = b2offset, j2n = b2offset * n; j2 != b2offset + b2size; ++j2, j2n += n) {
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____b + j0n + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                            data[m3___m2___k____b + j2n + j3] += den[_____m1___k____b + j1n + j0] * integral2;
                          }
                        }
                      }
                    }
                  } else {
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell
                      for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____b + j0n + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                          }
                        }
                      }
                    }
                  }

                }
              }
            }
          }
        }

      }
    }
  }

  delete[] diskdata;
  RefMatrix out(new PMatrix1e(coulomb_real_space->ft()));
  return out;
}


RefMatrix PMP2::coulomb_runtime() const {
  const double gamma = geom_->gamma();
  RefMatrix density(new PMatrix1e(coeff_->form_density_rhf()));
  const complex<double>* den = density->data()->front();

  RefMatrix coulomb_real_space(new PMatrix1e(union_geom_));
  complex<double>* data = coulomb_real_space->data()->front();

  const int s = geom_->S();
  const int l = geom_->L();

  /// oooo loop starts here.
  {
    const size_t allocsize_o = eri_obs_->max_num_int();
    double* diskdata_o = new double[allocsize_o];

    long file_position = 0l;
    size_t mcnt = 0lu;
    for (int m1 = -s; m1 <= s; ++m1) {
      for (int m2 = 0; m2 <= l; ++m2) { // use bra-ket symmetry!!!
        for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++mcnt) {

          const int k = geom_->K();
          const size_t b = density->blocksize();
          const size_t qb = coulomb_real_space->blocksize();

          const int n = geom_->nbasis();
          const int q = union_geom_->nbasis();

          // for density.
          const int m2___m3___k____b = (m2 - m3 + k) * b;
          const int _____m1___k____b = (   - m1 + k) * b;

          {
            const int m1________k____qb = (m1      + k) * qb;
            const int m3___m2___k____qb = (m3 - m2 + k) * qb;

            eri_obs_->get_block(file_position, eri_obs_->num_int_each(mcnt), diskdata_o);
            file_position += eri_obs_->num_int_each(mcnt);
            const double* cdata = diskdata_o;

            const int size = eri_obs_->basissize(); // number of shells
            for (int i0 = 0; i0 != size; ++i0) {
              const int b0offset = eri_obs_->offset(i0);
              const int b0size = eri_obs_->nbasis(i0);

              for (int i1 = 0; i1 != size; ++i1) {
                const int b1offset = eri_obs_->offset(i1);
                const int b1size = eri_obs_->nbasis(i1);

                for (int i2 = 0; i2 != size; ++i2) {
                  const int b2offset = eri_obs_->offset(i2);
                  const int b2size = eri_obs_->nbasis(i2);

                  for (int i3 = 0; i3 != size; ++i3) {
                    const int b3offset = eri_obs_->offset(i3);
                    const int b3size = eri_obs_->nbasis(i3);

                    const double integral_bound = eri_obs_->schwarz(((m1      + k) * size + i0) * size + i1)
                                                * eri_obs_->schwarz(((m3 - m2 + k) * size + i2) * size + i3);
                    const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                    if (skip_schwarz) continue;

                    if (m2 != 0) {
                      for (int j0 = b0offset, j0q = b0offset * q; j0 != b0offset + b0size; ++j0, j0q += q) { // center unit cell
                        for (int j1 = b1offset, j1n = b1offset * n; j1 != b1offset + b1size; ++j1, j1n += n) {
                          for (int j2 = b2offset, j2q = b2offset * q; j2 != b2offset + b2size; ++j2, j2q += q) {
                            for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                              const double integral2 = *cdata + *cdata;
                              data[m1________k____qb + j0q + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                              data[m3___m2___k____qb + j2q + j3] += den[_____m1___k____b + j1n + j0] * integral2;
                            }
                          }
                        }
                      }
                    } else {
                      for (int j0 = b0offset, j0q = b0offset * q; j0 != b0offset + b0size; ++j0, j0q += q) { // center unit cell
                        for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {
                          for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {
                            for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                              const double integral2 = *cdata + *cdata;
                              data[m1________k____qb + j0q + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                            }
                          }
                        }
                      }
                    }

                  }
                }
              }
            }
          }

        }
      }
    }
    delete[] diskdata_o;
  }
  /// oooo loop ends.

  /// ooco loop starts here.
  {
    const size_t allocsize_cs = eri_cabs_->max_num_int();
    double* diskdata_cs = new double[allocsize_cs];

    long file_position = 0l;
    size_t mcnt = 0lu;
    for (int m1 = -s; m1 <= s; ++m1) {
      for (int m2 = -l; m2 <= l; ++m2) { // *NOT* using bra-ket symmetry!!!
        for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++mcnt) {

          const int k = geom_->K();
          const size_t b = density->blocksize();
          const size_t qb = coulomb_real_space->blocksize();

          const int n = geom_->nbasis();
          const int q = union_geom_->nbasis();

          // for density.
          const int m2___m3___k____b = (m2 - m3 + k) * b;
          const int _____m1___k____b = (   - m1 + k) * b;

          {
            const int m1________k____qb = (m1      + k) * qb;
            const int m3___m2___k____qb = (m3 - m2 + k) * qb;

            eri_cabs_->get_block(file_position, eri_cabs_->num_int_each(mcnt), diskdata_cs);
            file_position += eri_cabs_->num_int_each(mcnt);
            const double* cdata = diskdata_cs;

            const int size_i = eri_cabs_->size_i();
            for (int i0 = 0; i0 != size_i; ++i0) { // OBS
              const int b0offset = eri_cabs_->offset_i(i0);
              const int b0size = eri_cabs_->basis_i(i0)->nbasis();

              const int size_a = eri_cabs_->size_a();
              for (int i1 = 0; i1 != size_a; ++i1) { // CABS!!!
                const int b1offset = eri_cabs_->offset_a(i1) + geom_->nbasis();
                const int b1size = eri_cabs_->basis_a(i1)->nbasis();

                const int size_j = eri_cabs_->size_j();
                for (int i2 = 0; i2 != size_j; ++i2) { // OBS
                  const int b2offset = eri_cabs_->offset_j(i2);
                  const int b2size = eri_cabs_->basis_j(i2)->nbasis();

                  const int size_b = eri_cabs_->size_b();
                  for (int i3 = 0; i3 != size_b; ++i3) { // OBS
                    const int b3offset = eri_cabs_->offset_b(i3);
                    const int b3size = eri_cabs_->basis_b(i3)->nbasis();

                    const double integral_bound = eri_cabs_->schwarz_ia(((m1      + k) * size_i + i0) * size_a + i1)
                                                * eri_cabs_->schwarz_jb(((m3 - m2 + k) * size_j + i2) * size_b + i3);
                    const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                    if (skip_schwarz) continue;

                    for (int j0 = b0offset, j0q = b0offset * q; j0 != b0offset + b0size; ++j0, j0q += q) { // center unit cell
                      for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____qb + j0q + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                          }
                        }
                      }
                    }

                  }
                }
              }
            }
          }

        }
      }
    }
    delete[] diskdata_cs;
  }
  /// ooco loop ends.

  /// cooo loop starts here.
  {
    shared_ptr<PCompCABSFile<ERIBatch> >
      eri_cabs_d(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, false, false, false, "ERI CABS(j)"));
    eri_cabs_d->store_integrals();
    eri_cabs_d->reopen_with_inout();
    const size_t allocsize_cs = eri_cabs_d->max_num_int();
    double* diskdata_cs = new double[allocsize_cs];

    long file_position = 0l;
    size_t mcnt = 0lu;
    for (int m1 = -s; m1 <= s; ++m1) {
      for (int m2 = -l; m2 <= l; ++m2) { // *NOT* using bra-ket symmetry!!!
        for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++mcnt) {

          const int k = geom_->K();
          const size_t b = density->blocksize();
          const size_t qb = coulomb_real_space->blocksize();

          const int n = geom_->nbasis();
          const int q = union_geom_->nbasis();

          // for density.
          const int m2___m3___k____b = (m2 - m3 + k) * b;
          const int _____m1___k____b = (   - m1 + k) * b;

          {
            const int m1________k____qb = (m1      + k) * qb;
            const int m3___m2___k____qb = (m3 - m2 + k) * qb;

            eri_cabs_d->get_block(file_position, eri_cabs_d->num_int_each(mcnt), diskdata_cs);
            file_position += eri_cabs_d->num_int_each(mcnt);
            const double* cdata = diskdata_cs;

            const int size_i = eri_cabs_d->size_i();
            for (int i0 = 0; i0 != size_i; ++i0) { // CABS!!!!
              const int b0offset = eri_cabs_d->offset_i(i0) + geom_->nbasis();
              const int b0size = eri_cabs_d->basis_i(i0)->nbasis();

              const int size_a = eri_cabs_d->size_a();
              for (int i1 = 0; i1 != size_a; ++i1) { // OBS
                const int b1offset = eri_cabs_d->offset_a(i1);
                const int b1size = eri_cabs_d->basis_a(i1)->nbasis();

                const int size_j = eri_cabs_d->size_j();
                for (int i2 = 0; i2 != size_j; ++i2) { // OBS
                  const int b2offset = eri_cabs_d->offset_j(i2);
                  const int b2size = eri_cabs_d->basis_j(i2)->nbasis();

                  const int size_b = eri_cabs_d->size_b();
                  for (int i3 = 0; i3 != size_b; ++i3) { // OBS
                    const int b3offset = eri_cabs_d->offset_b(i3);
                    const int b3size = eri_cabs_d->basis_b(i3)->nbasis();

                    const double integral_bound = eri_cabs_d->schwarz_ia(((m1      + k) * size_i + i0) * size_a + i1)
                                                * eri_cabs_d->schwarz_jb(((m3 - m2 + k) * size_j + i2) * size_b + i3);
                    const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                    if (skip_schwarz) continue;

                    for (int j0 = b0offset, j0q = b0offset * q; j0 != b0offset + b0size; ++j0, j0q += q) { // center unit cell
                      for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____qb + j0q + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                          }
                        }
                      }
                    }

                  }
                }
              }
            }
          }

        }
      }
    }
    delete[] diskdata_cs;
  }
  /// ooco loop ends.

  /// cooo loop starts here.
  {
    shared_ptr<PCompCABSFile<ERIBatch> >
      eri_cabs_t(new PCompCABSFile<ERIBatch>(geom_, gamma, true, false, true, false, false, "ERI CABS(ja)"));
    eri_cabs_t->store_integrals();
    eri_cabs_t->reopen_with_inout();
    const size_t allocsize_cs = eri_cabs_t->max_num_int();
    double* diskdata_cs = new double[allocsize_cs];

    long file_position = 0l;
    size_t mcnt = 0lu;
    for (int m1 = -s; m1 <= s; ++m1) {
      for (int m2 = -l; m2 <= l; ++m2) { // *NOT* using bra-ket symmetry!!!
        for (int m3 = m2 - s; m3 <= m2 + s; ++m3, ++mcnt) {

          const int k = geom_->K();
          const size_t b = density->blocksize();
          const size_t qb = coulomb_real_space->blocksize();

          const int n = geom_->nbasis();
          const int q = union_geom_->nbasis();

          // for density.
          const int m2___m3___k____b = (m2 - m3 + k) * b;
          const int _____m1___k____b = (   - m1 + k) * b;

          {
            const int m1________k____qb = (m1      + k) * qb;
            const int m3___m2___k____qb = (m3 - m2 + k) * qb;

            eri_cabs_t->get_block(file_position, eri_cabs_t->num_int_each(mcnt), diskdata_cs);
            file_position += eri_cabs_t->num_int_each(mcnt);
            const double* cdata = diskdata_cs;

            const int size_i = eri_cabs_t->size_i();
            for (int i0 = 0; i0 != size_i; ++i0) { // CABS!!!!
              const int b0offset = eri_cabs_t->offset_i(i0) + geom_->nbasis();
              const int b0size = eri_cabs_t->basis_i(i0)->nbasis();

              const int size_a = eri_cabs_t->size_a();
              for (int i1 = 0; i1 != size_a; ++i1) { // CABS!!!
                const int b1offset = eri_cabs_t->offset_a(i1) + geom_->nbasis();
                const int b1size = eri_cabs_t->basis_a(i1)->nbasis();

                const int size_j = eri_cabs_t->size_j();
                for (int i2 = 0; i2 != size_j; ++i2) { // OBS
                  const int b2offset = eri_cabs_t->offset_j(i2);
                  const int b2size = eri_cabs_t->basis_j(i2)->nbasis();

                  const int size_b = eri_cabs_t->size_b();
                  for (int i3 = 0; i3 != size_b; ++i3) { // OBS
                    const int b3offset = eri_cabs_t->offset_b(i3);
                    const int b3size = eri_cabs_t->basis_b(i3)->nbasis();

                    const double integral_bound = eri_cabs_t->schwarz_ia(((m1      + k) * size_i + i0) * size_a + i1)
                                                * eri_cabs_t->schwarz_jb(((m3 - m2 + k) * size_j + i2) * size_b + i3);
                    const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                    if (skip_schwarz) continue;

                    for (int j0 = b0offset, j0q = b0offset * q; j0 != b0offset + b0size; ++j0, j0q += q) { // center unit cell
                      for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____qb + j0q + j1] += den[m2___m3___k____b + j3n + j2] * integral2;
                          }
                        }
                      }
                    }

                  }
                }
              }
            }
          }

        }
      }
    }
    delete[] diskdata_cs;
  }

  return coulomb_real_space;
}



