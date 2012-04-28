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
       const std::vector<std::shared_ptr<Atom> >& aux_atoms,  const std::vector<std::vector<int> >& aux_offsets, const double thr,
       const bool j2, const std::shared_ptr<DensityFit> df = std::shared_ptr<DensityFit>())
     : DensityFit(nbas, naux) {
       common_init(atoms, offsets, atoms, offsets, aux_atoms, aux_offsets, thr, j2);

       if (!j2) {
         if (!df) throw std::logic_error("ERI fit should be called with j2 = true, IF df is not provided!"); 
         std::unique_ptr<double[]> d(new double[naux_*naux_]); 
         std::copy(df->data_2index(), df->data_2index()+naux_*naux_, d.get());
         data2_ = move(d);
       }

    };
    ~ERIFit() {};
    
};

#endif

