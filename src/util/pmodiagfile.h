//
// Author : Toru Shiozaki
// Date   : September 2009
//

#pragma once
#ifndef __src_util_pmodiagfile_h
#define __src_util_pmodiagfile_h

#include <src/util/pfile.h>
#include <src/pscf/f77.h>

template<class T>
class PMODiagFile : public PFile<T> {
  protected:
    const int istart_;
    const int jstart_;
    const int ifence_;
    const int jfence_;

  public:
    PMODiagFile(const long, const int, const int, const int, const int, const int, const bool late_init = false); 
    ~PMODiagFile(); 

};


template<class T>
PMODiagFile<T>::PMODiagFile(const long fsize, const int k,
                            const int istrt, const int ifen,
                            const int jstrt, const int jfen,
                            const bool late_init)
 : PFile<T>(fsize, k, late_init), istart_(istrt), ifence_(ifen), jstart_(jstrt), jfence_(jfen) {


};


template<class T>
PMODiagFile<T>::~PMODiagFile() {
};


#endif

