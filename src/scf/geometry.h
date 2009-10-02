//
// Author: Toru Shiozaki
// Date  : May 2009
//

#ifndef __src_scf_geometry_h
#define __src_scf_geometry_h

#include <string>
#include <vector>
#include <src/scf/atom.h>
#include <src/scf/petite.h>
#include <boost/shared_ptr.hpp>

class Geometry {
  protected:
    bool spherical_;

    std::string input_;
    std::vector<boost::shared_ptr<Atom> > atoms_;
    std::vector<boost::shared_ptr<Atom> > cabs_atoms_;
    double nuclear_repulsion_;

    int nbasis_;
    int nocc_;
    int nfrc_;
    int ncabs_;
    int lmax_;
    int cabs_lmax_;
    std::vector<std::vector<int> > offsets_;
    std::vector<std::vector<int> > cabs_offsets_;

    const double compute_nuclear_repulsion();

    int level_;
    std::string basisfile_;
    std::string cabsfile_;

    std::string symmetry_;
    boost::shared_ptr<Petite> plist_;
    int nirrep_;

  public:
    Geometry(const std::string, const int level);
    ~Geometry();

    std::vector<boost::shared_ptr<Atom> > atoms() const { return atoms_; };
    std::vector<boost::shared_ptr<Atom> > cabs_atoms() const { return cabs_atoms_; };
    boost::shared_ptr<Atom> atoms(const unsigned int i) const { return atoms_[i]; };

    const int natom() const { return atoms_.size(); };
    const int nbasis() const { return nbasis_; };
    const int nocc() const { return nocc_; };
    const int nfrc() const { return nfrc_; };
    const int ncabs() const { return ncabs_; };
    const int lmax() const { return lmax_; };
    const int cabs_lmax() const { return cabs_lmax_; };
    const bool spherical() const { return spherical_; };
    const int nirrep() const { return nirrep_; };
    const std::string symmetry() const { return symmetry_; };
    virtual const double nuclear_repulsion() const { return nuclear_repulsion_; };
    const int level() const { return level_; };
    const std::string basisfile() const { return basisfile_; };
    const std::string cabsfile() const { return cabsfile_; };

    const std::vector<std::vector<int> > offsets() const { return offsets_; };
    const std::vector<std::vector<int> > cabs_offsets() const { return cabs_offsets_; };
    const std::vector<int> offset(const unsigned int i) const { return offsets_.at(i); };

    void print_atoms() const;
    boost::shared_ptr<Petite> plist() { return plist_; }; 
};

#endif

