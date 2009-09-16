//
// Author : Toru Shiozaki
// Date   : June 2009
//

#include <src/scf/symmat.h>
#include <algorithm>
#include <iostream>

using namespace std;

typedef boost::shared_ptr<Geometry> RefGeometry;
typedef boost::shared_ptr<Atom> RefAtom;
typedef boost::shared_ptr<Shell> RefShell;
typedef boost::shared_ptr<Petite> RefPetite;
typedef boost::shared_ptr<SymRotAbel> RefSymRotAbel;

SymMat::SymMat(const RefGeometry gm, const int iop) : Matrix1e(gm), petite_(gm->plist()) {

  RefSymRotAbel symrt(new SymRotAbel(petite_->symop(iop), gm->lmax(), gm->spherical()));
  symrot_ = symrt;

  const vector<RefAtom> atoms = gm->atoms(); 
  const int natom = atoms.size();
  vector<RefAtom>::const_iterator aiter, biter;

  vector<int> offsets(natom, 0);
  for (int i = 1; i != natom; ++i) offsets[i] = offsets[i - 1] + atoms[i - 1]->nbasis(); 

  const bool spherical = gm->spherical();

  for (int i = 0; i != natom; ++i) {
    const int tatom_num = petite_->sym_atommap(i)[iop];
    const RefAtom catom = atoms[i]; 
    const RefAtom tatom = atoms[tatom_num];

    const vector<RefShell> shells = catom->shells(); 
    const int nshell = shells.size();

    int coffset = offsets[i];
    int toffset = offsets[tatom_num];
    for (int j = 0; j != nshell; ++j) {
      const int ang = shells[j]->angular_number(); 
      const int size = spherical ? (2 * ang + 1) : ((ang + 1) * (ang + 2) / 2);
      const int nfunc = shells[j]->num_contracted();
      const vector<double> block = symrot_->primrot(ang);
      assert(block.size() == size * size);

      for (int k = 0; k != nfunc; ++k) {
        int index = 0;
        for (int ic = coffset; ic != coffset + size; ++ic) {
          for (int it = toffset; it != toffset + size; ++it, ++index) 
            data_[ic * nbasis_ + it] = block[index]; 
        }
        coffset += size; 
        toffset += size; 
      } 
    }
  }

}


SymMat::~SymMat() {

} 
