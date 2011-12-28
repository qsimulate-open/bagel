//
// Author: Toru Shiozaki
// Date  : April 2009
//
#define PITWOHALF 17.493418327624862

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <src/slater/slaterbatch.h>
#include <src/util/f77.h>
#include <src/rysint/macros.h>
#include <src/rysint/carsphlist.h>

using namespace std;

void SlaterBatch::compute() {
  bool swapped = false;
  const double zero = 0.0;
  const int zeroint = 0;
  const int unit = 1;

  bkup_ = new double[size_alloc_];
  if (yukawa_) 
    bkup2_ = new double[size_alloc_];
  
  const int size = size_alloc_;
  dcopy_(&size, &zero, &zeroint, data_, &unit);
  if (yukawa_) 
    dcopy_(&size, &zero, &zeroint, data2_, &unit);

  // perform VRR
  // data_ and data2_ will contain the intermediates: prim01{ prim23{ xyz{ } } } 
  if (yukawa_) {
    switch (rank_) {
      case 1: perform_USVRR1(); break;
      case 2: perform_USVRR2(); break;
#if 0
// I cannot believe my code
// See comments in svrr_template.cc 
      case 3: perform_USVRR3(); break;
      case 4: perform_USVRR4(); break;
      case 5: perform_USVRR5(); break;
      case 6: perform_USVRR6(); break;
      case 7: perform_USVRR7(); break;
      case 8: perform_USVRR8(); break;
      case 9: perform_USVRR9(); break;
      case 10: perform_USVRR10(); break;
      case 11: perform_USVRR11(); break;
      case 12: perform_USVRR12(); break;
      case 13: perform_USVRR13(); break;  
#endif
      default: perform_USVRR(); break;
    } 
  } else {
    switch (rank_) {
      case 1: perform_SVRR1(); break;
      case 2: perform_SVRR2(); break;
#if 0
// I cannot believe my code
// See comments in svrr_template.cc 
      case 3: perform_SVRR3(); break;
      case 4: perform_SVRR4(); break;
      case 5: perform_SVRR5(); break;
      case 6: perform_SVRR6(); break;
      case 7: perform_SVRR7(); break;
      case 8: perform_SVRR8(); break;
      case 9: perform_SVRR9(); break;  
      case 10: perform_SVRR10(); break;  
      case 11: perform_SVRR11(); break;  
      case 12: perform_SVRR12(); break;  
      case 13: perform_SVRR13(); break;
#endif
      default: perform_SVRR(); break;
    } 
  }

  // contract indices 01 
  // data will be stored in bkup_: cont01{ prim23{ xyz{ } } }
  {
    const int m = prim2size_ * prim3size_ * asize_ * csize_; 
    const int bkupsize = m * cont0size_ * cont1size_; 
    perform_contraction_new_outer(m, data_, prim0size_, prim1size_, bkup_, 
            basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_, 
            basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
    if (yukawa_)
      perform_contraction_new_outer(m, data2_, prim0size_, prim1size_, bkup2_, 
              basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_, 
              basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
  }
  // contract indices 23 
  // data will be stored in data_: cont01{ cont23{ xyz{ } } }
  {
    const int n = cont0size_ * cont1size_;
    const int datasize = n * cont2size_ * cont3size_ * asize_ * csize_; 
    perform_contraction_new_inner(n, bkup_, prim2size_, prim3size_, data_, 
            basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_, 
            basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
    if (yukawa_)
      perform_contraction_new_inner(n, bkup2_, prim2size_, prim3size_, data2_, 
              basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_, 
              basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
  }

  // HRR to indices 01
  // data will be stored in bkup_: cont01{ cont23{ xyzf{ xyzab{ } } } }
  {
    if (basisinfo_[1]->angular_number() != 0) { 
      const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr_.hrrfunc_call(hrr_index, contsize_ * csize_, data_, AB_, bkup_);
      if (yukawa_)
        hrr_.hrrfunc_call(hrr_index, contsize_ * csize_, data2_, AB_, bkup2_);
    } else {
      swapped = true;
    }
  }

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();

  int a = (ang0 + 1) * (ang0 + 2) / 2;
  int b = (ang1 + 1) * (ang1 + 2) / 2;
  int c = (ang2 + 1) * (ang2 + 2) / 2;
  int d = (ang3 + 1) * (ang3 + 2) / 2;
  const int asph = 2 * ang0 + 1;
  const int bsph = 2 * ang1 + 1;
  const int csph = 2 * ang2 + 1;
  const int dsph = 2 * ang3 + 1;


  // Cartesian to spherical 01 if necesarry
  // data will be stored in data_
  struct CarSphList carsphlist;
  if (spherical_) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_ * csize_;
    if (!swapped) {
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_); 
      if (yukawa_) 
        carsphlist.carsphfunc_call(carsphindex, nloops, bkup2_, data2_); 
    } else {
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
      if (yukawa_) 
        carsphlist.carsphfunc_call(carsphindex, nloops, data2_, bkup2_); 
    }
    swapped = (swapped ^ true); 
  }


  // Transform batch for the second HRR step
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if cartesian
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if spherical
  {
    const int m = spherical_ ? (asph * bsph) : (a * b);
    const int n = cont2size_ * cont3size_ * csize_;
    const int nloop = cont0size_ * cont1size_;
    int offset = 0;
    if (swapped) {
      for (int i = 0; i != nloop; ++i, offset += m * n) {
        mytranspose_(&data_[offset], &m, &n, &bkup_[offset]);
        if (yukawa_) 
          mytranspose_(&data2_[offset], &m, &n, &bkup2_[offset]);
      }
    } else {
      for (int i = 0; i != nloop; ++i, offset += m * n) {
        mytranspose_(&bkup_[offset], &m, &n, &data_[offset]);
        if (yukawa_) 
          mytranspose_(&bkup2_[offset], &m, &n, &data2_[offset]);
      }
    }
  }

  // HRR to indices 23
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if cartesian
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if spherical
  {
    if (basisinfo_[3]->angular_number() != 0) { 
      const int hrr_index = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
      if (swapped && spherical_)       hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, bkup_, CD_, data_);
      else if (swapped)                hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, bkup_, CD_, data_);
      else if (!swapped && spherical_) hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, data_, CD_, bkup_);
      else                             hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, data_, CD_, bkup_);

      if (yukawa_) { 
        if (swapped && spherical_)       hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, bkup2_, CD_, data2_);
        else if (swapped)                hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, bkup2_, CD_, data2_);
        else if (!swapped && spherical_) hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, data2_, CD_, bkup2_);
        else                             hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, data2_, CD_, bkup2_);
      }
    } else {
      swapped = (swapped ^ true); 
    }
  }

  // Cartesian to spherical 23 if necesarry
  // data will be stored in bkup_
  if (spherical_) {
    const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    const int nloops = contsize_ * asph * bsph;
    if (swapped) {
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, data2_, bkup2_); 
    } else {
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_); 
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, bkup2_, data2_); 
    }
    swapped = (swapped ^ true);
    a = asph;
    b = bsph;
    c = csph;
    d = dsph;
  }

  double *data_now = swapped ? bkup_ : data_;
  double *bkup_now = swapped ? data_ : bkup_;
  double *data_now_2 = swapped ? bkup2_ : data2_;
  double *bkup_now_2 = swapped ? data2_ : bkup2_;

  // Sort cont23 and xyzcd
  // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
  {
    const int nloop = a * b * cont0size_ * cont1size_;
    const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
    sort_.sortfunc_call(index, data_now, bkup_now, cont3size_, cont2size_, nloop, swap23_);
    if (yukawa_)
      sort_.sortfunc_call(index, data_now_2, bkup_now_2, cont3size_, cont2size_, nloop, swap23_);
  }

  // transpose batch
  // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } } 
  {
    const int m = c * d * cont2size_ * cont3size_;
    const int n = a * b * cont0size_ * cont1size_; 
    mytranspose_(data_now, &m, &n, bkup_now);
    if (yukawa_)
      mytranspose_(data_now_2, &m, &n, bkup_now_2);
  }

  // Sort cont01 and xyzab
  // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
  {
    const int nloop = c * d * cont2size_ * cont3size_;
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(index, data_now, bkup_now, cont1size_, cont0size_, nloop, swap01_);
    if (yukawa_)
      sort_.sortfunc_call(index, data_now_2, bkup_now_2, cont1size_, cont0size_, nloop, swap01_);
  }
  
  if (swapped) {
    dcopy_(&size, bkup_, &unit, data_, &unit); 
    if (yukawa_)
      dcopy_(&size, bkup2_, &unit, data2_, &unit); 
  }

  delete[] bkup_;
  if (yukawa_) delete[] bkup2_;
}


