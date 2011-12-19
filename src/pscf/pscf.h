//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pscf_h
#define __src_pscf_pscf_h

#include <memory>
#include <src/util/pcompfile.h>
#include <src/pscf/pgeometry.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/phcore.h>
#include <src/pscf/pfock.h>
#include <src/pscf/ptildex.h>
#include <src/util/pdata.h> 
#include <src/pscf/pcoeff.h>
#include <src/rysint/eribatch.h>

class PSCF {
  protected:
    const std::shared_ptr<PGeometry> geom_;
    const std::shared_ptr<PHcore> hcore_;
    const std::shared_ptr<POverlap> overlap_;
    std::shared_ptr<PTildeX> tildex_;
    std::shared_ptr<PMatrix1e> aodensity_;
    std::shared_ptr<PCoeff> coeff_;
    std::vector<double> schwarz_;

    virtual void init_schwarz();
    bool direct_;

    // created at constructor, deleted at destructor.
    double* eig_;
    void print_eig(const double*) const;
    std::shared_ptr<PCompFile<ERIBatch> > ao_eri_;

    const double obtain_energy(const PMatrix1e&, const PMatrix1e&, const PMatrix1e&);

  public:
    PSCF(const std::shared_ptr<PGeometry>);
    PSCF(const std::shared_ptr<PGeometry>, std::shared_ptr<PMatrix1e>);
    ~PSCF();

    virtual void compute();
    const std::shared_ptr<PGeometry> geom() { return geom_;} ;
    const std::shared_ptr<PCoeff> coeff() { return coeff_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    // Returning a copy of eig_. Returns std::vector<double>
    const std::vector<double> eig() const {
      const size_t length = (2*geom_->K()+1) * geom_->nbasis();
      const std::vector<double> out(eig_, eig_+length);
      return out;
    }
    std::shared_ptr<PCompFile<ERIBatch> > ao_eri() { return ao_eri_; };
};

#endif

