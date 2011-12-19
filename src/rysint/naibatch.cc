//
// Author: Toru Shiozaki
// Date  : April 2009
//
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <src/rysint/inline.h>
#include <src/rysint/naibatch.h>
#include <src/rysint/f77.h>
#include <src/rysint/macros.h>
#include <src/rysint/erirootlist.h>
#define PI 3.1415926535897932
#define SQRTPI2 0.886226925452758013649083741671

using namespace std;

typedef std::shared_ptr<Geometry> RefGeometry;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

NAIBatch::NAIBatch(const vector<RefShell> _info, const RefGeometry gm, const int L, const double A)
 :  RysInt(_info), geom_(gm), L_(L), A_(A) {

  const double integral_thresh = PRIM_SCREEN_THRESH; 
//const double integral_thresh = 0.0;

  // natom_ 
  natom_ = geom_->natom() * (2 * L + 1);

  // cartesian to spherical tranformation has not yet been implemented.
  const bool spherical = false;

  // swap 01 indices when needed: Larger angular momentum function comes first
  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number()) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();

  rank_ = ceil(0.5 * (ang0 + ang1 + 1));

  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);

  prim0size_ = basisinfo_[0]->num_primitive();
  prim1size_ = basisinfo_[1]->num_primitive();
  primsize_ = prim0size_ * prim1size_;
  cont0size_ = basisinfo_[0]->num_contracted();
  cont1size_ = basisinfo_[1]->num_contracted();
  contsize_ = cont0size_ * cont1size_;

  amax_ = ang0 + ang1;
  amin_ = max(ang0, ang1);
  amax1_ = amax_ + 1;

  asize_ = 0; 
  for (int i = amin_; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;

  asize_intermediate_ = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4; 
  asize_final_ = spherical_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_intermediate_;

  int cnt = 0;
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

  const unsigned int size_start = asize_ * primsize_; 
  size_final_ = asize_final_ * contsize_;
  const unsigned int size_intermediate = asize_intermediate_ * primsize_;
  size_alloc_ = max(max(size_start, size_intermediate), size_final_);
  data_ = new double[size_alloc_];
  fill(data_, data_ + size_alloc_, 0.0);
  assert(size_final_ <= size_alloc_);

  buff_ = new double[(rank_ * 2 + 6) * primsize_ * natom_];
  double* pointer = buff_; 
  p_ = pointer;     pointer += primsize_ * natom_ * 3;
  xp_ = pointer;    pointer += primsize_ * natom_;
  coeff_ = pointer; pointer += primsize_ * natom_;
  T_ = pointer;     pointer += primsize_ * natom_;

  screening_ = new int [primsize_ * natom_];
  screening_size_ = 0;

  vector<double>::const_iterator expi0, expi1, expi2, expi3;
  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;
  vector<RefAtom> atoms;
  // for periodic calculations
  if (L == 0) {
    atoms =  geom_->atoms();
  } else {
    vector<RefAtom> tmp =  geom_->atoms();
    for (int i = -L; i <= L; ++i) {
      vector<double> disp(3);
      disp[0] = disp[1] = 0.0;
      disp[2] = A * i; 
      for (int j = 0; j != tmp.size(); ++j) {
        RefAtom newatom(new Atom(*(tmp[j]), disp));
        atoms.push_back(newatom);
      }
    }
  }
  const double onepi2 = 1.0 / (PI * PI);
  const double sqrtpi = ::sqrt(PI);
  for (expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) { 
    for (expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1) { 
      for (vector<RefAtom>::const_iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++index) {
        double Z = static_cast<double>((*aiter)->atom_number()); 
        const double cxp = *expi0 + *expi1;
        xp_[index] = cxp;
        const double ab = *expi0 * *expi1; 
        const double cxp_inv = 1.0 / cxp;
        p_[index * 3    ] = (basisinfo_[0]->position(0) * *expi0 + basisinfo_[1]->position(0) * *expi1) * cxp_inv;
        p_[index * 3 + 1] = (basisinfo_[0]->position(1) * *expi0 + basisinfo_[1]->position(1) * *expi1) * cxp_inv;
        p_[index * 3 + 2] = (basisinfo_[0]->position(2) * *expi0 + basisinfo_[1]->position(2) * *expi1) * cxp_inv;
        const double Eab = ::exp(-(AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]) * (ab * cxp_inv) );
        coeff_[index] = - 2 * Z * PI * cxp_inv * Eab; 
        const double PCx = p_[index * 3    ] - (*aiter)->position(0);
        const double PCy = p_[index * 3 + 1] - (*aiter)->position(1);
        const double PCz = p_[index * 3 + 2] - (*aiter)->position(2);
        const double T = cxp * (PCx * PCx + PCy * PCy + PCz * PCz); 
        const double sqrtt = ::sqrt(T);
        const double ss = - coeff_[index] * ::pow(4.0 * ab * onepi2, 0.75) * (T > 1.0e-15 ? sqrtpi * ::erf(sqrtt) / sqrtt * 0.5 : 1.0); 
        if (ss > integral_thresh) {
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
        } else {
          T_[index] = -1.0;
          coeff_[index] = 0.0;
        }
      }
    }
  }

  roots_ = pointer; pointer += rank_ * primsize_ * natom_; 
  weights_ = pointer;
  fill(weights_, weights_ + rank_ * primsize_ * natom_, 0.0);

  int ps = (int)primsize_ * natom_; 

  if(rank_ != 1) {
    rysroot_(T_, roots_, weights_, &rank_, &ps);
  } else {
    for (int i = 0; i != ps; ++i) {
      const double t = T_[i];
      const double sqrtt = ::sqrt(t);
      weights_[i] = t < 1.0e-10 ? 1.0 : inline_erf(sqrtt) * SQRTPI2 / sqrtt;
      roots_[i] = t < 1.0e-10 ? 0.5 : (weights_[i] - ::exp(-t)) / t * 0.5 / weights_[i];
    }
  }

}


NAIBatch::~NAIBatch() {
  delete[] data_;
  delete[] buff_;
  delete[] screening_;
}

