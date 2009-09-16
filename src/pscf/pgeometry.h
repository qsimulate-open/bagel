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

    double pnuclear_repulsion_;

  public:
    PGeometry(const std::string, const int);
    ~PGeometry() {};

    const int L() const { return L_; };
    const int S() const { return S_; };
    const int K() const { return K_; };
    const double A() const { return A_; };

    const double compute_pnuclear_repulsion();

    const double nuclear_repulsion() const { return pnuclear_repulsion_; }; 
};

#endif

