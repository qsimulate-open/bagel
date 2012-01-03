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
#include <memory>
#include <src/df/df.h>
#include <src/util/input.h>

class Geometry {
  protected:
    // Spherical or Cartesian basis set.
    bool spherical_;

    // Input file name.
    std::string input_;

    // Atoms, which contains basis-set info also.
    std::vector<std::shared_ptr<Atom> > atoms_;
    std::vector<std::shared_ptr<Atom> > aux_atoms_;
    bool aux_merged_;

    // Nuclear repulsion energy.
    double nuclear_repulsion_;
    // Computes the nuclear repulsion energy.
    const double compute_nuclear_repulsion();

    // Some shared info for basis sets.
    int nbasis_;
    int nocc_;
    int nfrc_;
    int naux_;
    int lmax_;
    int aux_lmax_;
    std::vector<std::vector<int> > offsets_;
    std::vector<std::vector<int> > aux_offsets_;

    int level_;
    std::string basisfile_;
    std::string auxfile_;

    // Symmetry can be used for molecular calculation.
    std::string symmetry_;
    std::shared_ptr<Petite> plist_;
    int nirrep_;

    // integral screening
    double schwarz_thresh_;

    // for DF calculations
    std::shared_ptr<DensityFit> df_;

    // for R12 calculations
    double gamma_;

  public:
    Geometry(const std::string, const int level);
    Geometry(const std::shared_ptr<InputData> inpt);
    ~Geometry();

    // Returns shared pointers of Atom objects, which contains basis-set info.
    std::vector<std::shared_ptr<Atom> > atoms() const { return atoms_; };
    std::vector<std::shared_ptr<Atom> > aux_atoms() const { return aux_atoms_; };
    std::shared_ptr<Atom> atoms(const unsigned int i) const { return atoms_[i]; };

    // Returns a constant
    const int natom() const { return atoms_.size(); };
    const size_t nbasis() const { return nbasis_; };
    const size_t nocc() const { return nocc_; };
    const size_t nfrc() const { return nfrc_; };
    const size_t naux() const { return naux_; };
    const int lmax() const { return lmax_; };
    const int aux_lmax() const { return aux_lmax_; };
    const bool spherical() const { return spherical_; };
    const int nirrep() const { return nirrep_; };
    const double gamma() const {return gamma_; };
    const std::string symmetry() const { return symmetry_; };
    virtual const double nuclear_repulsion() const { return nuclear_repulsion_; };
    const int level() const { return level_; };
    const std::string basisfile() const { return basisfile_; };
    const std::string auxfile() const { return auxfile_; };
    const double schwarz_thresh() const { return schwarz_thresh_; };

    // TODO for some reasons needed now in CASSCF
    void set_nocc(const int i) { nocc_ = i; };
    void set_basis(const int i) { nbasis_ = i; };
    void set_ncore(const int i) { nfrc_ = i; };
    int num_count_ncore(); // also set nfrc_
    int num_count_full_valence_nocc();

    // The position of the specific funciton in the basis set.
    const std::vector<std::vector<int> > offsets() const { return offsets_; };
    const std::vector<std::vector<int> > aux_offsets() const { return aux_offsets_; };
    const std::vector<int> offset(const unsigned int i) const { return offsets_.at(i); };
    const std::vector<int> aux_offset(const unsigned int i) const { return aux_offsets_.at(i); };

    // Printing out some info
    void print_atoms() const;

    // Returns the Petite list.
    std::shared_ptr<Petite> plist() { return plist_; }; 

    // Rerurns DF data
    std::shared_ptr<DensityFit> df() { return df_; };

    // In R12 methods, we need to construct a union of OBS and CABS.
    // Currently, this is done by creating another object and merge OBS and CABS into atoms_.
    // After this, compute_nuclear_repulsion() should not be called.
    // Not undo-able.
    void merge_obs_aux() {
      aux_merged_ = true;
      atoms_.insert(atoms_.end(), aux_atoms_.begin(), aux_atoms_.end());
      for (std::vector<std::vector<int> >::iterator iter = aux_offsets_.begin(); iter != aux_offsets_.end(); ++iter) {
        for (std::vector<int>::iterator citer = iter->begin(); citer != iter->end(); ++citer) {
          *citer += nbasis_;
        }
      }
      offsets_.insert(offsets_.end(), aux_offsets_.begin(), aux_offsets_.end());
      nbasis_ += naux_;
    };
};

#endif

