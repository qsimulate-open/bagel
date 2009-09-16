//
// Author : Toru Shiozaki
// Date   : June 2009
//

#ifndef __src_scf_symrot_h
#define __src_scf_symrot_h

#include <vector>

class SymRotAbel {
  protected:
    std::vector<std::vector<double> > primrot_; 

  public:
    std::vector<double> primrot(const int i) const { return primrot_[i]; };

    SymRotAbel(const std::vector<double>&, const int, const bool);
    ~SymRotAbel();


};

#endif

