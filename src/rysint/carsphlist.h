//
// Author: Toru Shiozaki
// Date  : May 2009
//
  
#ifndef __rysint_carsphlist_h
#define __rysint_carsphlist_h

#include <src/rysint/macros.h>
 
struct CarSphList {
  CarSphList();
  ~CarSphList();
     
  static void carsph_00(const int, const double*, double*);
  static void carsph_10(const int, const double*, double*);
  static void carsph_20(const int, const double*, double*);
  static void carsph_30(const int, const double*, double*);
  static void carsph_40(const int, const double*, double*);
  static void carsph_50(const int, const double*, double*);
  static void carsph_60(const int, const double*, double*);
  static void carsph_11(const int, const double*, double*);
  static void carsph_21(const int, const double*, double*);
  static void carsph_31(const int, const double*, double*);
  static void carsph_41(const int, const double*, double*);
  static void carsph_51(const int, const double*, double*);
  static void carsph_61(const int, const double*, double*);
  static void carsph_22(const int, const double*, double*);
  static void carsph_32(const int, const double*, double*);
  static void carsph_42(const int, const double*, double*);
  static void carsph_52(const int, const double*, double*);
  static void carsph_62(const int, const double*, double*);
  static void carsph_33(const int, const double*, double*);
  static void carsph_43(const int, const double*, double*);
  static void carsph_53(const int, const double*, double*);
  static void carsph_63(const int, const double*, double*);
  static void carsph_44(const int, const double*, double*);
  static void carsph_54(const int, const double*, double*);
  static void carsph_64(const int, const double*, double*);
  static void carsph_55(const int, const double*, double*);
  static void carsph_65(const int, const double*, double*);
  static void carsph_66(const int, const double*, double*);

  void carsphfunc_call(const unsigned int i, const int a0, const double* a1, double* a2) {
    return (carsphfunc[i])(a0, a1, a2);
  };

  void (*carsphfunc[ANG_HRR_END * ANG_HRR_END])(const int, const double*, double*);
};

#endif
