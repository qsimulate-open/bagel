//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __NEWINT_SRC_CASSCF_QVEC
#define __NEWINT_SRC_CASSCF_QVEC

#include <memory>
#include <src/df/df.h>
#include <src/scf/coeff.h>
#include <src/fci/fci.h> // 2RDM and half-transformed integrals
#include <src/casscf/rotfile.h>

class Qvec : public QFile {
  protected:

  public:
    Qvec(const int n, const int m, std::shared_ptr<DensityFit> df, std::shared_ptr<Coeff> c, const size_t nclosed,
         std::shared_ptr<FCI> fci);
    Qvec(const QFile& a) : QFile(a) {};
    ~Qvec() {};

};


#endif
