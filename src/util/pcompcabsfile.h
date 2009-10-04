//
// Author : Toru Shiozaki
// Date   : October 2009
//

#pragma once
#ifndef __src_util_pcompcabsfile_h
#define __src_util_pcompcabsfile_h

#include <src/util/pcompfile.h>

template<class T>
class PCompCABSFile : public PCompFile<T> {

  protected:

  public:
    PCompCABSFile(boost::shared_ptr<PGeometry>, const bool late_init = false, const std::string jobname = "source");


};


template<class T>
PCompCABSFile<T>::PCompCABSFile(boost::shared_ptr<PGeometry> pg, const bool late_init, const std::string jobname)
 : PCompFile<T>(pg, late_init, jobname) {

};

#endif

