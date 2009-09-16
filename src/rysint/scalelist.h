//
// Author: Toru Shiozaki
// Date  : May 2009
// 

#ifndef __src_rysint_scalelist_h
#define __src_rysint_scalelist_h

#include <src/rysint/macros.h>

struct ScaleList {
  public:
    ScaleList();
    ~ScaleList();

    static void scale_data_1(double*, const double*, const double, const double*, const int);
    static void scale_data_2(double*, const double*, const double, const double*, const int);
    static void scale_data_3(double*, const double*, const double, const double*, const int);
    static void scale_data_4(double*, const double*, const double, const double*, const int);
    static void scale_data_5(double*, const double*, const double, const double*, const int);
    static void scale_data_6(double*, const double*, const double, const double*, const int);
    static void scale_data_7(double*, const double*, const double, const double*, const int);
    static void scale_data_8(double*, const double*, const double, const double*, const int);
    static void scale_data_9(double*, const double*, const double, const double*, const int);
    static void scale_data_10(double*, const double*, const double, const double*, const int);
    static void scale_data_11(double*, const double*, const double, const double*, const int);
    static void scale_data_12(double*, const double*, const double, const double*, const int);
    static void scale_data_13(double*, const double*, const double, const double*, const int);

    void (*scalefunc[RYS_MAX + 1])(double*, const double*, const double, const double*, const int);
};

#endif

