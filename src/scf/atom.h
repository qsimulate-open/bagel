//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __scf_atom_h
#define __scf_atom_h

#include <vector>
#include <string>
#include <src/scf/shell.h>
#include <src/scf/atommap.h>
#include <memory>

class Atom {
  protected:
    bool spherical_;

    std::string name_;
    std::vector<double> position_;
    std::vector<std::shared_ptr<Shell> > shells_; 
    int atom_number_;
    int nbasis_;
    int lmax_;

    AtomMap atommap_;
     
  public:
    Atom(const bool, const std::string, const std::vector<double>&, const std::string);
    Atom(const Atom&, const std::vector<double>&);
    Atom(const Atom&, const double*);
    ~Atom();

    const std::string name() const { return name_; };
    const int atom_number() const { return atom_number_;};
    const std::vector<double> position() const { return position_; };
    const double position(const unsigned int i) const { return position_[i]; };
    const std::vector<std::shared_ptr<Shell> > shells() const { return shells_; };
    const int nshell() const { return shells_.size(); };

    const int nbasis() const { return nbasis_; };
    const int lmax() const { return lmax_; };
    const bool spherical() const { return spherical_; };

    void print_basis() const;
    void print() const;
};

#endif

