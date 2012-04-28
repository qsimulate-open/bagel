//
// Author : Toru Shiozaki
// Date   : April 2012
//

#ifndef __SRC_DF_FIT_H
#define __SRC_DF_FIT_H

#include <src/df/df.h>

class ERIFit : public DensityFit {
  protected:
    void compute_batch(std::unique_ptr<double[]>& data2_, std::vector<std::shared_ptr<Shell> >& input,
                       const int j0st, const int j0fen, const int j1st, const int j1fen, const int naux) {
      ERIBatch eribatch(input, 0.0);
      eribatch.compute();
      const double* eridata = eribatch.data();
      for (int j0 = j0st; j0 != j0fen; ++j0)
        for (int j1 = j1st; j1 != j1fen; ++j1, ++eridata)
          data2_[j1+j0*naux] = data2_[j0+j1*naux] = *eridata; 
    };
  public:
    ERIFit(const int nbas, const int naux,
       const std::vector<std::shared_ptr<Atom> >& atoms,  const std::vector<std::vector<int> >& offsets,
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr)
     : DensityFit(nbas, naux) {
       common_init(atoms, offsets, aux_atoms, aux_offsets, thr);
    };
    ~ERIFit() {};
    
};

#endif

