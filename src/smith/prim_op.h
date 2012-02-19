//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __SRC_SMITH_PRIM_OP_H
#define __SRC_SMITH_PRIM_OP_H

#include <memory>

namespace SMITH {

static void sort_indices4(const std::unique_ptr<double[]>& unsorted, std::unique_ptr<double[]>& sorted,
                          const int a,const int b,const int c,const int d, // according to unsorted
                          const int i,const int j,const int k,const int l,
                          const double factor = 1.0) {
  int id[4];
  int jd[4] = {a, b, c, d};

  long iall=0;
  for(int j0=0;j0<a;++j0){
    id[0]=j0;
    for(int j1=0;j1<b;++j1){
      id[1]=j1;
      for(int j2=0;j2<c;++j2){
        id[2]=j2;
        for (int j3=0;j3<d;++j3,++iall){
          id[3]=j3;
          long ib=id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]));
          sorted[ib]=unsorted[iall]*factor;
        }
      }
    } 
  } 
};

}

#endif
