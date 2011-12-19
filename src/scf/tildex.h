//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_tildex_h
#define __src_scf_tildex_h

#include <src/scf/overlap.h>
#include <src/scf/matrix1e.h>
#include <memory>

class TildeX : public Matrix1e {
  protected:

  public:
    TildeX(const std::shared_ptr<Overlap>);
    ~TildeX();

};

#endif
