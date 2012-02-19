//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_PRIM_OP_H
#define __SRC_SMITH_PRIM_OP_H

#include <memory>

namespace SMITH {

static void sort_indices2(const std::unique_ptr<double[]>& unsorted, std::unique_ptr<double[]>& sorted,
                          const int b, const int a, // according to unsorted
                          const int i, const int j, const double factor = 1.0) {
  if (i==0) {
    const int j0max = a*b;
    for (int j0=0; j0<a*b; ++j0) sorted[j0]=unsorted[j0]*factor;
  } else { 
    int id[2];
    int jd[2] = {b, a};
    long iall=0;
    for(int j0=0;j0<(int)a;++j0){
      id[1]=j0;
      for(int j1=0;j1<(int)b;++j1,++iall){
        id[0]=j1;
        long ib=id[i]+jd[i]*id[j];
        sorted[ib]=unsorted[iall]*factor;
      }
    }
  }
} 


// CAUTION :: I have changed the convention from that in mpqc.
static void sort_indices4(const std::unique_ptr<double[]>& unsorted, std::unique_ptr<double[]>& sorted,
                          const int d,const int c,const int b,const int a, // according to unsorted
                          const int i,const int j,const int k,const int l,
                          const double factor = 1.0) {
  int id[4];
  int jd[4] = {d, c, b, a};

  long iall=0;
  for(int j0=0;j0<a;++j0){
    id[3]=j0;
    for(int j1=0;j1<b;++j1){
      id[2]=j1;
      for(int j2=0;j2<c;++j2){
        id[1]=j2;
        for (int j3=0;j3<d;++j3,++iall){
          id[0]=j3;
          long ib=id[i]+jd[i]*(id[j]+jd[j]*(id[k]+jd[k]*id[l]));
          sorted[ib]=unsorted[iall]*factor;
        }
      }
    } 
  } 
};

}

#endif
