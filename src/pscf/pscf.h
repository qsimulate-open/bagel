//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pscf_h
#define __src_pscf_pscf_h

#include <boost/shared_ptr.hpp>
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
    const boost::shared_ptr<PGeometry> geom_;
    const boost::shared_ptr<PHcore> hcore_;
    const boost::shared_ptr<POverlap> overlap_;
    boost::shared_ptr<PTildeX> tildex_;
    boost::shared_ptr<PMatrix1e> aodensity_;
    boost::shared_ptr<PCoeff> coeff_;
    std::vector<double> schwarz_;

    virtual void init_schwarz();
    bool direct_;

    // created at constructor, deleted at destructor.
    double* eig_;
    void print_eig(const double*) const;
    boost::shared_ptr<PCompFile<ERIBatch> > ao_eri_;

    const double obtain_energy(const PMatrix1e&, const PMatrix1e&, const PMatrix1e&);

  public:
    PSCF(const boost::shared_ptr<PGeometry>);
    ~PSCF();

    virtual void compute();
    const boost::shared_ptr<PGeometry> geom() { return geom_;} ;
    const boost::shared_ptr<PCoeff> coeff() { return coeff_; };
    const std::vector<double>& schwarz() const { return schwarz_; };

    // Returning a copy of eig_. Returns std::vector<double>
    const std::vector<double> eig() const {
      const size_t length = (2*geom_->K()+1) * geom_->nbasis();
      const std::vector<double> out(eig_, eig_+length);
      return out;
    }
    boost::shared_ptr<PCompFile<ERIBatch> > ao_eri() { return ao_eri_; };
};

#endif

