//
// Author : Toru Shiozaki
// Date   : June 2009
//

#ifndef __src_scf_symmat_h
#define __src_scf_symmat_h

#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <src/scf/petite.h>
#include <src/scf/symrot.h>

class SymMat : public Matrix1e {
  protected:
    std::shared_ptr<SymRotAbel> symrot_; 
    std::shared_ptr<Petite> petite_; 
    void computebatch(const std::vector<std::shared_ptr<Shell> >&, const int, const int, const int) {};

  public:
    SymMat(const std::shared_ptr<Geometry>, const int);
    ~SymMat(); 

};

#endif
