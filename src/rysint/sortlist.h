//
// Author: Toru Shiozaki
// Date  : April 2009
//
//
// stand alone class for HRR functions.
// can be called from outside through sortfunc_call.
//
  
#ifndef __src_rysint_sortlist_h
#define __src_rysint_sortlist_h

#include <src/rysint/macros.h>
 
struct SortList {
  SortList(const bool);
  ~SortList();
     
  static void sort_indices_00(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_01(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_02(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_03(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_04(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_05(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_06(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_11(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_12(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_13(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_14(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_15(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_16(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_22(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_23(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_24(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_25(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_26(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_33(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_34(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_35(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_36(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_44(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_45(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_46(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_55(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_56(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_66(double*, const double*, const int, const int, const int, const bool);

  static void sort_indices_00_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_01_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_02_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_03_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_04_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_05_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_06_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_11_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_12_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_13_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_14_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_15_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_16_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_22_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_23_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_24_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_25_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_26_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_33_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_34_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_35_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_36_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_44_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_45_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_46_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_55_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_56_sph(double*, const double*, const int, const int, const int, const bool);
  static void sort_indices_66_sph(double*, const double*, const int, const int, const int, const bool);

  void sortfunc_call(const unsigned int i, double* a1, const double* a2, const int a3, const int a4, const int a5, const bool a6) {
    return (sortfunc[i])(a1, a2, a3, a4, a5, a6);
  };

  void (*sortfunc[ANG_HRR_END * ANG_HRR_END])(double*, const double*, const int, const int, const int, const bool);
};

#endif
