//header




#ifndef __BAGEL_ZFCI_ZPHASE_H
#define __BAGEL_ZFCI_ZPHASE_H
#include <src/zfci/zknowles.h>

namespace bagel {

std::complex<double>* ZKnowlesHandy::phase_factor(std::shared_ptr<const MOFile> jop, int nelec, int ij) const {

    std::complex<double> complex_moe[ij*ij];
//generating complex matrix from mofile real matrix
  if(nelec==2) {
  //temporary fix to allow for non-complex mo2e.
    std::cout << std::endl << std::endl << "printing complex_m02e" << std::endl << std::endl;
    for (int j=0;j<(ij);j++) {
      std::cout << j << std::setw(1);
      for (int i=0; i<(ij); i++) {
        complex_moe[i+j*ij] = std::complex<double>(jop->mo2e(i,j), 0.0);
        std::cout << std::setprecision(1) << complex_moe[i+j*ij] << std::setw(1);
      }
      std::cout << std::endl;
    }
  }
  else if(nelec==1) {
  //temporary fix to allow for non-complex mo1e.
    std::cout << std::endl << std::endl << "printing complex_m01e" << std::endl << std::endl;
    for (int j=0;j<(ij);j++) {
      std::cout << j << std::setw(1);
      for (int i=0; i<(ij); i++) {
        complex_moe[i] = std::complex<double>(jop->mo1e(i), 0.0);
        std::cout << std::setprecision(1) << complex_moe[i+j*ij] << std::setw(1);
      }
      std::cout << std::endl;
    }
  }
//introducing phase factor to 2e integrals to check for bugs in complex portion of zfci
//random numbers for use in phase factor
  std::complex<double> random[ij];
  for (int j=0;j<ij;j++) {
    double a = rand() % 10 + 1;
    random[j] = std::complex<double>(0,a);
  }
//declaring and filling intermediate w/ 0.0
  std::complex<double> intermediate[ij*ij];
  for (int j=0; j<(ij*ij); j++) intermediate[j] =0.0;

//generating diagonal phase factor matrix
  std::complex<double> phase_factor[ij*ij];
  std::cout << std::endl << std::endl << "printing phase factor" << std::endl << std::endl;
  for (int j=0;j<(ij);j++) {
    std::cout << j << std::setw(1);
    for (int i=0; i<(ij); i++) {
      if (i==j) phase_factor[i+j*ij] = random[i];
      else phase_factor[i+j*ij] = 0.0;
      assert(phase_factor[i+j*ij].real()<1e-8);
      std::cout << std::setprecision(1) << phase_factor[i+j*ij] << std::setw(1);
    }
    std::cout << std::endl;
  }
// (A*U)
  zgemm3m_("n", "n", ij, ij, ij, 1.0, complex_moe, ij, phase_factor, ij, 0.0, intermediate, ij);

  std::cout << std::endl << std::endl << "printing intermediate" << std::endl << std::endl;
  for (int j=0; j<(ij); j++) {
    std::cout << j << std::setw(1);
    for (int i=0; i<(ij); i++) {
      std::cout << std::setprecision(1) << intermediate[i+j*ij] << std::setw(1);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

// Udagger*(A*U)
  zgemm3m_("c", "n", ij, ij, ij, 1.0, phase_factor, ij, intermediate, ij, 0.0, complex_moe, ij);

//end of phase factor
//there has to be a better way to do this, but this is temporary anyway
  std::unique_ptr<std::complex<double>[]> out;
  out = std::unique_ptr<std::complex<double>[]> (new std::complex<double>[ij*ij]);
  for (int i=0;i<(ij*ij);i++) out[i] = complex_moe[i];
  return out.get();
}
}
#endif
