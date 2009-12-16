//
// Author: Toru Shiozaki
// Date  : April 2009
//
#define PITWOHALF 17.493418327624862
#define PIMHALF 0.564189583547756
#define SQRTPI2 0.886226925452758013649083741671

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <src/slater/slaterbatch.h>
#include <src/slater/srootlist.h>
#include <src/slater/f77.h>
#include <src/rysint/macros.h>
#include <src/rysint/inline.h>
#include <src/slater/sinline.h>

using namespace std;
using namespace boost;

#define MIN_EXPONENT -200.0
#define GM1_THRESH 1.0e-15
#define GM1_THRESH_MAX 1.0e+15

typedef boost::shared_ptr<Shell> RefShell;

SlaterBatch::SlaterBatch(const vector<RefShell> _info, const double max_density, const double gmm, const bool yukawa)
:  RysInt(_info), gamma_(gmm), yukawa_(yukawa) {

  const double integral_thresh = std::max(1.0e-30, PRIM_SCREEN_THRESH * 1.0e-10);// / max_density;

  // swap 01 indices when needed: Larger angular momentum function comes first
  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number()) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }
  // swap 23 indices when needed
  if (basisinfo_[2]->angular_number() < basisinfo_[3]->angular_number()) {
    swap(basisinfo_[2], basisinfo_[3]);
    swap23_ = true;
  } else {
    swap23_ = false;
  }

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();

  rank_ = ceil(0.5 * (ang0 + ang1 + ang2 + ang3)) + 1; // for slater!!

  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);

  AB_[0] = ax - bx; 
  AB_[1] = ay - by; 
  AB_[2] = az - bz; 
  CD_[0] = cx - dx; 
  CD_[1] = cy - dy; 
  CD_[2] = cz - dz; 

  prim0size_ = basisinfo_[0]->num_primitive();
  prim1size_ = basisinfo_[1]->num_primitive();
  prim2size_ = basisinfo_[2]->num_primitive();
  prim3size_ = basisinfo_[3]->num_primitive();
  primsize_ = prim0size_ * prim1size_ * prim2size_ * prim3size_;
  cont0size_ = basisinfo_[0]->num_contracted();
  cont1size_ = basisinfo_[1]->num_contracted();
  cont2size_ = basisinfo_[2]->num_contracted();
  cont3size_ = basisinfo_[3]->num_contracted();
  contsize_ = cont0size_ * cont1size_ * cont2size_ * cont3size_;

  amax_ = ang0 + ang1;
  cmax_ = ang2 + ang3;
  amin_ = ang0;
  cmin_ = ang2;
  amax1_ = amax_ + 1;
  cmax1_ = cmax_ + 1;

  asize_ = 0; 
  for (int i = amin_; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;
  csize_ = 0; 
  for (int i = cmin_; i <= cmax_; ++i) csize_ += (i + 1) * (i + 2) / 2;

  const int asize_final = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4;
  const int csize_final = (ang2 + 1) * (ang2 + 2) * (ang3 + 1) * (ang3 + 2) / 4;

  const int asize_final_sph = spherical_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_final;
  const int csize_final_sph = spherical_ ? (2 * ang2 + 1) * (2 * ang3 + 1) : csize_final;

  int cnt = 0;
  for (int i = cmin_; i <= cmax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) { 
      for (int iy = 0; iy <= i - iz; ++iy) { 
        const int ix = i - iy - iz;
        if (ix >= 0) {
          cmapping_[ix + cmax1_ * (iy + cmax1_ * iz)] = cnt;
          ++cnt;
        }
      }
    }
  }
  cnt = 0;
  for (int j = amin_; j <= amax_; ++j) {
    for (int jz = 0; jz <= j; ++jz) { 
      for (int jy = 0; jy <= j - jz; ++jy) { 
        const int jx = j - jy - jz;
        if (jx >= 0){
          amapping_[jx + amax1_ * (jy + amax1_ * jz)] = cnt; 
          ++cnt;
        }
      }
    }
  }
  
  const unsigned int size_start = asize_ * csize_ * primsize_; 
  const unsigned int size_intermediate = asize_final * csize_ * contsize_;
  const unsigned int size_intermediate2 = asize_final_sph * csize_final * contsize_;
  size_final_ = asize_final_sph * csize_final_sph * contsize_;
  size_alloc_ = max(size_start, max(size_intermediate, size_intermediate2));
  data_ = new double[size_alloc_];
  data2_ = new double[size_alloc_];

  buff_ = new double[(rank_ * 2 + 12) * primsize_];
  double* pointer = buff_; 
  p_ = pointer;     pointer += primsize_ * 3;
  q_ = pointer;     pointer += primsize_ * 3;
  xp_ = pointer;    pointer += primsize_;
  xq_ = pointer;    pointer += primsize_;
  coeff_ = pointer; pointer += primsize_;
  coeffy_ = pointer;pointer += primsize_;
  T_ = pointer;     pointer += primsize_;
  U_ = pointer;     pointer += primsize_;

  vector<double>::const_iterator expi0, expi1, expi2, expi3;
  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();
  const vector<double> exp2 = basisinfo_[2]->exponents();
  const vector<double> exp3 = basisinfo_[3]->exponents();

  vector<double> Ecd_save(prim2size_ * prim3size_); 
  vector<double> qx_save(prim2size_ * prim3size_);
  vector<double> qy_save(prim2size_ * prim3size_);
  vector<double> qz_save(prim2size_ * prim3size_);

  const double minexp0 = *min_element(exp0.begin(), exp0.end());
  const double minexp1 = *min_element(exp1.begin(), exp1.end());
  const double minexp2 = *min_element(exp2.begin(), exp2.end());
  const double minexp3 = *min_element(exp3.begin(), exp3.end());
  const double min_ab = minexp0 * minexp1;
  const double min_cd = minexp2 * minexp3;

  // minimum distance between two lines (AB and CD)
  const double x_ab_cd = AB_[1] * CD_[2] - AB_[2] * CD_[1]; 
  const double y_ab_cd = AB_[2] * CD_[0] - AB_[0] * CD_[2]; 
  const double z_ab_cd = AB_[0] * CD_[1] - AB_[1] * CD_[0]; 
  const double x_n_ac = x_ab_cd * (ax - cx);
  const double y_n_ac = y_ab_cd * (ay - cy);
  const double z_n_ac = z_ab_cd * (az - cz);
  const double innerproduct = x_n_ac + y_n_ac + z_n_ac; 
  const double norm_ab_cd_sq = x_ab_cd * x_ab_cd + y_ab_cd * y_ab_cd + z_ab_cd * z_ab_cd;
  const double min_pq_sq = norm_ab_cd_sq == 0.0 ? 0.0 : innerproduct * innerproduct / norm_ab_cd_sq; 

  indexpair23_.reserve(prim2size_ * prim3size_);

  const double r01_sq = AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]; 
  const double r23_sq = CD_[0] * CD_[0] + CD_[1] * CD_[1] + CD_[2] * CD_[2];

  {
    const double cxp_min = minexp0 + minexp1; 
    const double cxp_inv_min = 1.0 / cxp_min; 
    const double min_Eab = ::exp(-r01_sq * min_ab * cxp_inv_min);
    int index23 = 0;
    for (expi2 = exp2.begin(); expi2 != exp2.end(); ++expi2) { 
      for (expi3 = exp3.begin(); expi3 != exp3.end(); ++expi3, ++index23) { 
        const double cxq = *expi2 + *expi3;
        const double cd = *expi2 * *expi3;
        const double cxq_inv = 1.0 / cxq;

        if (-r23_sq * (cd * cxq_inv) < MIN_EXPONENT) continue;
        Ecd_save[index23] = ::exp(-r23_sq * (cd * cxq_inv) );
        qx_save[index23] = (cx * *expi2 + dx * *expi3) * cxq_inv;
        qy_save[index23] = (cy * *expi2 + dy * *expi3) * cxq_inv;
        qz_save[index23] = (cz * *expi2 + dz * *expi3) * cxq_inv;

        indexpair23_.push_back(make_tuple(index23, *expi2, *expi3));
      }
    }
  }

  int index = 0;
  int index01 = 0;
  fill(coeff_, coeff_ + primsize_, 0.0);
  fill(coeffy_, coeffy_ + primsize_, 0.0);
  fill(T_, T_ + primsize_, -1.0);
  fill(U_, U_ + primsize_, 1.0e-100);
  screening_ = new int[primsize_];
  screening_size_ = 0;

  const double cxq_min = minexp2 + minexp3; 
  const double cxq_inv_min = 1.0 / cxq_min; 
  const double min_Ecd = ::exp(-r23_sq * min_cd * cxq_inv_min);
  const double twogamma = 2.0 / gamma_;
  for (expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) { 
    for (expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1, ++index01) { 
      const double cxp = *expi0 + *expi1;
      const double ab = *expi0 * *expi1; 
      const double cxp_inv = 1.0 / cxp;
      if (-r01_sq * (ab * cxp_inv) < MIN_EXPONENT) continue;
      const double Eab = ::exp(-r01_sq * (ab * cxp_inv) );
      const double coeff_half = 2 * Eab * PITWOHALF;
      const double px = (ax * *expi0 + bx * *expi1) * cxp_inv;
      const double py = (ay * *expi0 + by * *expi1) * cxp_inv;
      const double pz = (az * *expi0 + bz * *expi1) * cxp_inv;

      const int index_base = prim2size_ * prim3size_ * index01;

      for (vector<tuple<int, double, double> >::const_iterator expi23 = indexpair23_.begin(); 
                                                               expi23 != indexpair23_.end(); ++expi23) {
          const int index23 = expi23->get<0>();
          const int index = index_base + index23;
          const double exp2value = expi23->get<1>(); 
          const double exp3value = expi23->get<2>(); 
          const double cxq = exp2value + exp3value;
          xp_[index] = cxp;
          xq_[index] = cxq; 
          const double cxpxq = cxp * cxq;
          const double onepqp_q = 1.0 / (::sqrt(cxp + cxq) * cxpxq);
          coeffy_[index] = Ecd_save[index23] * coeff_half * onepqp_q;
          const double rho = cxpxq / (cxp + cxq);
          const double xpq = qx_save[index23] - px;
          const double ypq = qy_save[index23] - py;
          const double zpq = qz_save[index23] - pz;
          const double pq_sq = xpq * xpq + ypq * ypq + zpq * zpq;
          const double U = 0.25 * gamma_ * gamma_ / rho;
          if (U - gamma_ * sqrt(pq_sq) < MIN_EXPONENT) continue;
          if (coeffy_[index] * exp(U - gamma_ * sqrt(pq_sq)) < integral_thresh) continue;

          const double T = rho * pq_sq;
          coeff_[index] = coeffy_[index] * twogamma * U;
          const int index3 = index * 3;
          p_[index3] = px; 
          p_[index3 + 1] = py;
          p_[index3 + 2] = pz;
          q_[index3] =     qx_save[index23]; 
          q_[index3 + 1] = qy_save[index23];
          q_[index3 + 2] = qz_save[index23];
          T_[index] = T; 
          U_[index] = U; 
          screening_[screening_size_] = index;
          ++screening_size_;
      }
    }
  }

  roots_ = pointer; pointer += rank_ * primsize_; 
  weights_ = pointer;
  fill(weights_, weights_ + primsize_, 0.0);

  // determine the quadrature grid
  if (rank_ == -1) {
    const double prefac = SQRTPI2 * 0.5;
    for (int i = 0; i != screening_size_; ++i) {
      const int ii = screening_[i];
      const double ct = T_[ii];
      const double cu = U_[ii];
      const double sqct = ::sqrt(ct);
      const double sqcu = ::sqrt(cu);
      const double lambda = sqct + sqcu;
      const double kappa = sqcu - sqct;
      double gm1, g0;

      if (fabs(ct) < 1.0e-10) {
        const double experfc_lambda = experfc(lambda);
        gm1 = SQRTPI2 / sqcu * experfc_lambda;
        g0  = 1 - (cu + cu) * gm1;
      } else {
        if (kappa < 0.0) {
          const double experfc_lambda = experfc(lambda);
          const double erfc_kappa = experfcm(kappa);
          const double factor1 = ::exp(-ct) * prefac;
          const double factor2 = ::exp(kappa * kappa - ct) * prefac;
          gm1 = (factor2 * erfc_kappa + factor1 * experfc_lambda) / sqcu;
          g0 =  (factor2 * erfc_kappa - factor1 * experfc_lambda) / sqct;
        } else {
          const double experfc_lambda = experfc(lambda);
          const double experfc_kappa  = experfc(kappa);
          const double factor = ::exp(-ct) * prefac;
          gm1 = factor / sqcu * (experfc_kappa + experfc_lambda);
          g0 =  factor / sqct * (experfc_kappa - experfc_lambda);
        }
      }

      if (gm1 > GM1_THRESH) {
        weights_[ii] = gm1; 
        roots_[ii] = g0 / gm1;
      } else {
        weights_[ii] = 0.0; 
        roots_[ii] = 0.0;
      }
    }
  } else if (rank_ == -2) {
    const double prefac = SQRTPI2 * 0.5;
    for (int i = 0; i != screening_size_; ++i) {
      const int ii = screening_[i];
      const double ct = T_[ii];
      const double cu = U_[ii];
      const double sqct = ::sqrt(ct);
      const double sqcu = ::sqrt(cu);
      const double lambda = sqct + sqcu;
      const double kappa = sqcu - sqct;
      double gm1, g0, g1, g2;

      if (fabs(ct) < 1.0e-10) {
        const double experfc_lambda = experfc(lambda);
        const double cu2 = cu + cu;
        gm1 = SQRTPI2 / sqcu * experfc_lambda;
        g0  = 1.0 - cu2 * gm1;
        g1  = (1.0 - cu2 * g0) * 0.33333333333333333;
        g2  = (1.0 - cu2 * g1) * 0.2;
      } else {
        if (kappa < 0.0) {
          const double experfc_lambda = experfc(lambda);
          const double erfc_kappa = experfcm(kappa);
          const double expct = ::exp(-ct);
          const double factor1 = expct * prefac;
          const double factor2 = ::exp(kappa * kappa - ct) * prefac;
          gm1 = (factor2 * erfc_kappa + factor1 * experfc_lambda) / sqcu;
          g0 = (factor2 * erfc_kappa - factor1 * experfc_lambda) / sqct;
          const double cu2 = cu + cu;
          const double one2t = 0.5 / ct;
          g1 = (g0 + cu2 * gm1 - expct) * one2t;  
          g2 = (3.0 * g1 + cu2 * g0 - expct) * one2t;  
        } else {
          const double experfc_lambda = experfc(lambda);
          const double experfc_kappa  = experfc(kappa);
          const double expct = ::exp(-ct);
          const double factor = expct * prefac;
          gm1 = factor / sqcu * (experfc_kappa + experfc_lambda);
          g0 = factor / sqct * (experfc_kappa - experfc_lambda);
          const double cu2 = cu + cu;
          const double one2t = 0.5 / ct;
          g1 = (g0 + cu2 * gm1 - expct) * one2t;  
          g2 = (3.0 * g1 + cu2 * g0 - expct) * one2t;  
        }
      }

      if (gm1 < GM1_THRESH) {
        const int offset =  ii + ii;
        roots_[offset] = 0.0;
        roots_[offset + 1] = 0.0;
        weights_[offset] = 0.0;
        weights_[offset + 1] = 0.0;
      } else {
        double x[2];
        double w[1];
        x[0] = g0 / gm1;
        const double sigma11 = g1 - x[0] * g0; 
        const double sigma12 = g2 - x[0] * g1; 
        x[1] = sigma12 / sigma11 - x[0]; 
        w[0] = sigma11 / gm1;

        const double dd = x[0] + x[1];
        double g = x[0] - x[1];
        g = ::sqrt(g * g + 4.0 * w[0]);
        const double s = (dd - g) * 0.5;
        const double c = (dd + g) * 0.5;
        const double bb = w[0];
        const double p = x[0] - s;
        const double f = x[0] - c;
        const double lp1 = bb / (p * p + bb);
        const double lp2 = bb / (f * f + bb);

        const int offset =  ii + ii;
        roots_[offset] = s;
        roots_[offset + 1] = c;
        weights_[offset] = lp1 * gm1;
        weights_[offset + 1] = lp2 * gm1;
      }
    }
  } else {
    int ps = (int)primsize_; 
    struct SRootList root;
    root.srootfunc_call(rank_, T_, U_, roots_, weights_, &ps); 
  }


#if 0
  cout << endl;
  int jk = 0;
  for (int i = 0; i != primsize_; ++i) { 
    cout << T_[i] << " " << U_[i] << endl;
    for (int j = 0; j != rank_; ++j, ++jk) 
      cout << fixed << setprecision(17) << setw(25) << roots_[jk] << " " << weights_[jk] << endl;
    cout << endl;
  }
#endif
}


SlaterBatch::~SlaterBatch() {
  delete[] data2_;
  delete[] data_;
  delete[] buff_;
  delete[] screening_;

}

