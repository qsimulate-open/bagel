//
// Author : Toru Shiozaki
// Date   : July 2009
//

#ifndef __src_pscf_pgeometry_h
#define __src_pscf_pgeometry_h

#include <string>
#include <src/scf/geometry.h>

class PGeometry : public Geometry {
  protected:
    int L_; // Namur cutoff L
    int S_; // Namur cutoff S
    int K_; // number of K points (= number of unit cell in the *half* first Brillouin zone)  
    double A_; // fundamental vector

    // Computes the nuclear repulsion energy per unit cell.
    double pnuclear_repulsion_;
    const double compute_pnuclear_repulsion();

  public:
    PGeometry(const std::string, const int);
    ~PGeometry() {};

    // Some constants for periodic calculations.
    const int L() const { return L_; };
    const int S() const { return S_; };
    const int K() const { return K_; };
    const double A() const { return A_; };

    // Returns nuclear repulsion energies.
    const double nuclear_repulsion() const { return pnuclear_repulsion_; }; 
};

#endif

