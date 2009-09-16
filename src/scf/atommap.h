//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_atommap_h
#define __src_scf_atommap_h

#include <map>
#include <string>

struct AtomMap {
  public:
    AtomMap(); 
    ~AtomMap();

    std::map<std::string, int> atommap;
    std::map<std::string, int> angmap;
};

#endif
