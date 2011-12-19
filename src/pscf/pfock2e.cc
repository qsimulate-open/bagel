//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <complex>
#include <iostream>
#include <src/pscf/pfock.h>
#include <src/rysint/eribatch.h>
#include <src/macros.h>
#include <fstream>
#include <cassert>

typedef std::complex<double> Complex;
typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;


using namespace std;

void PFock::pfock_two_electron_part() {

  assert(!direct());

  const int shift = sizeof(int) * 4;

  vector<size_t> tmp = file_->num_int_each(); 
  size_t allocsize = *max_element(tmp.begin(), tmp.end()); 
  double* diskdata = new double[allocsize];

  long file_position = 0l;
  size_t mcnt = 0lu;
  for (int m1 = -S2_; m1 <= S2_; ++m1) {
    for (int m2 = 0; m2 <= L(); ++m2) { // use bra-ket symmetry!!!
      for (int m3 = m2 - S2_; m3 <= m2 + S2_; ++m3, ++mcnt) {

        const bool m2nonzero = m2 != 0;
        const bool m2m1in = abs(m2 - m1) <= K();
        const bool m3in = abs(m3) <= K();
        const bool mmmin = m2m1in && m3in;

        const int k = K();
        const size_t b = blocksize_;
        const int n = nbasis_;

        const int m1________k____b = (m1      + k) * b;
        const int m1___m2___k____b = (m1 - m2 + k) * b;
        const int m2___m1___k____b = (m2 - m1 + k) * b;
        const int m2___m3___k____b = (m2 - m3 + k) * b;
        const int m3___m2___k____b = (m3 - m2 + k) * b;
        const int m3________k____b = (m3      + k) * b;
        const int _____m1___k____b = (   - m1 + k) * b;
        const int _____m3___k____b = (   - m3 + k) * b;

        {
          file_->get_block(file_position, file_->num_int_each(mcnt), diskdata);
          file_position += file_->num_int_each(mcnt);
          const double* cdata = diskdata;
          complex<double>* density = density_->data()->front();
          complex<double>* data = data_->front();
          
          const int size = file_->basissize(); // number of shells
          for (int i0 = 0; i0 != size; ++i0) {
            const int b0offset = file_->offset(i0); 
            const int b0size = file_->nbasis(i0);

            for (int i1 = 0; i1 != size; ++i1) {
              const int b1offset = file_->offset(i1);
              const int b1size = file_->nbasis(i1);

              for (int i2 = 0; i2 != size; ++i2) {
                const int b2offset = file_->offset(i2);
                const int b2size = file_->nbasis(i2);

                for (int i3 = 0; i3 != size; ++i3) {
                  const int b3offset = file_->offset(i3);
                  const int b3size = file_->nbasis(i3);

                  const double integral_bound = file_->schwarz(((m1      + k) * size + i0) * size + i1)
                                              * file_->schwarz(((m3 - m2 + k) * size + i2) * size + i3);
                  const bool skip_schwarz = integral_bound < SCHWARZ_THRESH;
                  if (skip_schwarz) continue;

                  if (m2nonzero && mmmin){
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell 
                      for (int j1 = b1offset, j1n = b1offset * n; j1 != b1offset + b1size; ++j1, j1n += n) {  
                        for (int j2 = b2offset, j2n = b2offset * n; j2 != b2offset + b2size; ++j2, j2n += n) {  
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {  
                            const double integral1 = *cdata;
                            const double integral2 = integral1 + integral1;
                            data[m1________k____b + j0n + j1] += density[m2___m3___k____b + j3n + j2] * integral2;
                            data[m3___m2___k____b + j2n + j3] += density[_____m1___k____b + j1n + j0] * integral2;
                            data[m3________k____b + j0n + j3] -= density[m2___m1___k____b + j1n + j2] * integral1; 
                            data[m1___m2___k____b + j2n + j1] -= density[_____m3___k____b + j3n + j0] * integral1; 
                          }
                        }
                      }
                    }
                  } else if (m2nonzero) {
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell 
                      for (int j1 = b1offset, j1n = b1offset * n; j1 != b1offset + b1size; ++j1, j1n += n) {  
                        for (int j2 = b2offset, j2n = b2offset * n; j2 != b2offset + b2size; ++j2, j2n += n) {  
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {  
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____b + j0n + j1] += density[m2___m3___k____b + j3n + j2] * integral2;
                            data[m3___m2___k____b + j2n + j3] += density[_____m1___k____b + j1n + j0] * integral2;
                          }
                        }
                      }
                    }
                  } else if (mmmin) {
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell 
                      for (int j1 = b1offset, j1n = b1offset * n; j1 != b1offset + b1size; ++j1, j1n += n) {  
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {  
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {  
                            const double integral1 = *cdata;
                            const double integral2 = integral1 + integral1;
                            data[m1________k____b + j0n + j1] += density[m2___m3___k____b + j3n + j2] * integral2;
                            data[m3________k____b + j0n + j3] -= density[m2___m1___k____b + j1n + j2] * integral1; 
                          }
                        }
                      }
                    }
                  } else {
                    for (int j0 = b0offset, j0n = b0offset * n; j0 != b0offset + b0size; ++j0, j0n += n) { // center unit cell 
                      for (int j1 = b1offset                    ; j1 != b1offset + b1size; ++j1          ) {  
                        for (int j2 = b2offset                    ; j2 != b2offset + b2size; ++j2          ) {  
                          for (int j3 = b3offset, j3n = b3offset * n; j3 != b3offset + b3size; ++j3, j3n += n, ++cdata) {  
                            const double integral2 = *cdata + *cdata;
                            data[m1________k____b + j0n + j1] += density[m2___m3___k____b + j3n + j2] * integral2;
                          }
                        }
                      }
                    }
                  }

                }
              }
            }
          }
        }

      }
    }
  } 

  delete[] diskdata;

}  

