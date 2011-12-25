//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

#include <iomanip>
#include <stdexcept>
#include <src/fci/fci.h>
#include <src/rysint/eribatch.h>
#include <boost/algorithm/combination.hpp>

using namespace std;
using namespace boost;

typedef std::shared_ptr<Atom> RefAtom;
typedef std::shared_ptr<Shell> RefShell;

static const int unit = 1;

//
// generate initial vectors
//   - bits: bit patterns of low-energy determinants
//   - nspin: #alpha - #beta
//   - out:

void FCI::generate_guess(const int nspin, const int nstate, std::shared_ptr<Dvec> out) {

  // TODO currently this is only for singlet states, sorry...
  assert(nspin == 0);

  int ndet = num_state_*10;
  start_over:
  vector<pair<int, int> > bits = detseeds(ndet);

  if (nspin == 0) {
    // in this case, easy. The singlet combinations are made for open-shell singlet bits
    int oindex = 0;
    vector<int> done;
    for (auto iter = bits.begin(); iter != bits.end(); ++iter) {
      const int alpha = iter->first;
      const int beta = iter->second;
      const int open_bit = (alpha^beta);

      // check if this orbital configuration is already used
      if (find(done.begin(), done.end(), open_bit) != done.end()) continue;
      done.push_back(open_bit);

      const int common = (alpha & beta);
      const int nalpha = numofbits(alpha^common);

      assert(nspin == 0 && numofbits(beta^common) == nalpha);

      vector<int> open = bit_to_numbers(open_bit);
      int init_alpha = common;
      for (int i=0; i!=nalpha; ++i) init_alpha &= (1<<open[i]);

      // take a linear combination to make a vector singlet coupled.
      int icnt = 0;
      do {
        int ialpha = common;
        int ibeta = common;
        for (int i=0; i!=nalpha; ++i) ialpha ^= (1<<open[i]);
        for (int i=nalpha; i!=open.size(); ++i) ibeta  ^= (1<<open[i]);
#if 0
      cout << "     string" << setw(3) << oindex << ":   alpha " <<
            print_bit(ialpha) << " beta " << print_bit(ibeta) << endl;
#endif
        const double sign = pow(-1.0, numofbits((init_alpha^ialpha)/2));
        out->data(oindex)->element(lexical<1>(ibeta), lexical<0>(ialpha)) = sign;
        ++icnt;
      } while (next_combination(open.begin(), open.begin()+nalpha, open.end()));

      // scale to make the vector normalized
      const double factor = 1.0/sqrt(static_cast<double>(icnt));
      const int size = out->data(oindex)->size();
      dscal_(&size, &factor, out->data(oindex)->first(), &unit);

      cout << "     guess " << setw(3) << oindex << ":   closed " <<
            setw(20) << left << print_bit(common) << " open " << setw(20) << print_bit(open_bit) << right << endl;

      ++oindex;
      if (oindex == nstate) break;
    }
    if (oindex < nstate) {
      out->zero(); 
      ndet *= 4;
      goto start_over;
    }
//    throw runtime_error("increase the number of determinants in the FCI guess generation");
    cout << endl;
  }

}


//
// returns seed determinants for initial guess
//
vector<pair<int, int> > FCI::detseeds(const int ndet) {
  multimap<double, pair<int,int> > tmp;
  for (int i = 0; i != ndet; ++i) tmp.insert(make_pair(-1.0e10*(1+i), make_pair(0,0)));

  double* diter = denom_->first();
  for (auto aiter = stringa_.begin(); aiter != stringa_.end(); ++aiter) {
    for (auto biter = stringb_.begin(); biter != stringb_.end(); ++biter, ++diter) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        tmp.insert(make_pair(din, make_pair(*biter, *aiter)));
        tmp.erase(tmp.begin());
      } 
    }
  }
  assert(tmp.size() == ndet || ndet > stringa_.size()*stringb_.size());
  vector<pair<int, int> > out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter) {
    out.push_back(iter->second);
#if 0
    cout << print_bit(iter->second.first) << " " << print_bit(iter->second.second) << " " << setprecision(10) << iter->first << endl;
#endif
  }
  return out;
}

//
// averaged diagonal elements as defined in Knowles & Handy (1989) Compt. Phys. Comm. 
//
void FCI::const_denom() {

  vector<double> jop, kop, fk;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  fk.resize(norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e(j, j, i, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e(j, i, j, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    fk[i] = 0.0;
    for (int j = 0; j != norb_; ++j) {
      fk[i] += kop[i*norb_+j];
    }
  }
  shared_ptr<Civec> tmp(new Civec(stringb_.size(), stringa_.size()));
  denom_ = tmp;
  const int nspin = numofbits(stringa_.front()) - numofbits(stringb_.front());
  const int nspin2 = nspin*nspin;

  double* iter = denom_->first();
  for (auto ia = stringa_.begin(); ia != stringa_.end(); ++ia) {
    for (auto ib = stringb_.begin(); ib != stringb_.end(); ++ib, ++iter) {
      unsigned int iabit1 = *ia;
      unsigned int ibbit1 = *ib;
      const int nopen = numofbits(iabit1^ibbit1);
      const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
      *iter = 0.0;
      for (int i = 0; i != norb_; ++i, (iabit1 >>= 1), (ibbit1 >>= 1)) {
        const int nia = (iabit1&1);
        const int nib = (ibbit1&1);
        const int niab = nia + nib;
        const int Ni = (nia ^ nib);
        unsigned int iabit2 = *ia;
        unsigned int ibbit2 = *ib;
        for (int j = 0; j != i; ++j, (iabit2 >>= 1), (ibbit2 >>= 1)) {
          const int nja = (iabit2&1);
          const int njb = (ibbit2&1);
          const int Nj = (nja ^ njb);
          const int addj = niab * (nja + njb); 
          *iter += jop[j+norb_*i] * 2.0 * addj - kop[j+norb_*i] * (F*Ni*Nj + addj);
        }
        *iter += (jop_->mo1e(i,i) + fk[i]) * niab - kop[i+norb_*i] * 0.5 * (Ni - niab*niab);
      }
    }
  }
}


void FCI::print_header() const {
  cout << "  ---------------------------" << endl;
  cout << "        FCI calculation      " << endl;
  cout << "  ---------------------------" << endl << endl;
}


#if 0
// I don't care the efficiency at all!!
// >>>>>>>>>>>>  ALL INCORE!!! <<<<<<<<<<<
void FCI::create_Jiiii() {
  // first compute all the AO integrals in core

  const int nocc = norb();
  const int size = jop_->basis().size(); // number of shells
  const int nbasis = geom_->nbasis();
  const size_t aointsize = nbasis*nbasis*nbasis*nbasis; 
  double* aobuff = new double[aointsize];
  // TODO this can be thrown away.
  double* first = new double[nbasis*nbasis*nbasis*nocc];

  // some stuffs for blas
  const int unit = 1;
  const double one = 1.0;
  const double zero = 0.0;
  double* cdata = ref_->coeff()->data() + ncore()*nbasis;

  cout << "  - AO integrals are computed and stored in core" << endl << endl;

  // one electron part
  {
    shared_ptr<Fock> fock0(new Fock(geom_, ref_->hcore()));
    if (ncore() != 0) {
      shared_ptr<Matrix1e> den(new Matrix1e(ref_->coeff()->form_core_density_rhf()));
      shared_ptr<Fock> fock1(new Fock(geom_, fock0, den, ref_->shwarz()));
      const double core_energy = (*den * (*ref_->hcore()+*fock1)).trace();
      set_core_energy(core_energy);
      fock0 = fock1;
    }
    fock0->symmetrize();
    dgemm_("n","n",&nbasis,&nocc,&nbasis,&one,fock0->data(),&nbasis,cdata,&nbasis,&zero,aobuff,&nbasis);
  }
  jop_->mo1e().resize(nocc*nocc);
  dgemm_("t","n",&nocc,&nocc,&nbasis,&one,cdata,&nbasis,aobuff,&nbasis,&zero,jop_->mo1e_ptr(),&nocc);

  for (int i0 = 0; i0 != size; ++i0) {
    const int b0offset = jop_->offset(i0); 
    const int b0size = jop_->basis(i0)->nbasis();
    for (int i1 = 0; i1 != size; ++i1) {
      const int b1offset = jop_->offset(i1);
      const int b1size = jop_->basis(i1)->nbasis();
      for (int i2 = 0; i2 != size; ++i2) {
        const int b2offset = jop_->offset(i2);
        const int b2size = jop_->basis(i2)->nbasis();
        for (int i3 = 0; i3 != size; ++i3) {
          const int b3offset = jop_->offset(i3); 
          const int b3size = jop_->basis(i3)->nbasis();
          vector<RefShell> input;
          input.push_back(jop_->basis(i3));
          input.push_back(jop_->basis(i2));
          input.push_back(jop_->basis(i1));
          input.push_back(jop_->basis(i0));

          ERIBatch eribatch(input, 1.0);
          eribatch.compute();
          
          const double* eridata = eribatch.data();
          // what a bad code!
          int ioff = 0;
          for (int i = b0offset; i != b0offset+b0size; ++i) {
            for (int j = b1offset; j != b1offset+b1size; ++j) {
              for (int k = b2offset; k != b2offset+b2size; ++k) {
                for (int l = b3offset; l != b3offset+b3size; ++l, ++ioff) {
                  aobuff[l+nbasis*(k+nbasis*(j+nbasis*i))] = eridata[ioff];
                }
              }
            }
          } 
          
        }
      }
    }
  }

  const int n = aointsize;
//cout << ddot_(&n, aobuff, &unit, aobuff, &unit) << endl; 
  const int nn = nbasis*nbasis;
  const int nnn = nbasis*nbasis*nbasis;
  dgemm_("n","n",&nnn,&nocc,&nbasis, &one, aobuff,&nnn,cdata,&nbasis, &zero,first,&nnn);

  for (int i = 0; i != nocc; ++i) {
    dgemm_("n","n",&nn,&nocc,&nbasis, &one, first+nnn*i,&nn,cdata,&nbasis, &zero,aobuff+nn*nocc*i,&nn);
  }

  const int mm = nocc*nocc;
  mytranspose_(aobuff,&nn,&mm,first);

  const int nmm = nbasis * mm; 
  dgemm_("n","n",&nmm,&nocc,&nbasis,&one,first,&nmm,cdata,&nbasis,&zero,aobuff,&nmm);

  for (int i = 0; i != nocc; ++i) {
    dgemm_("n","n",&mm,&nocc,&nbasis,&one,aobuff+nmm*i,&mm,cdata,&nbasis,&zero,first+mm*nocc*i,&mm);
  }
  delete[] aobuff;

  // storing unpacked integrals
  jop_->mo1e_unpacked().resize(mm);
  copy(jop_->mo1e().begin(), jop_->mo1e().end(), jop_->mo1e_unpacked().begin());
  jop_->mo2e_unpacked().resize(mm*mm);
  copy(first, first+mm*mm, jop_->mo2e_unpacked_ptr());

  // mo2e is compressed
  const int sizeij = nocc*(nocc+1)/2;
  jop_->set_sizeij(sizeij);
  jop_->mo2e().resize(sizeij*sizeij);
  // TODO very bad
  int ijkl = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j <= i; ++j) {
      const int ijo = (j + i*nocc)*nocc*nocc;
      for (int k = 0; k != nocc; ++k) {
        for (int l = 0; l <= k; ++l, ++ijkl) {
          *jop_->mo2e_ptr(ijkl) = first[l+k*nocc+ijo]; 
        }
      }
    }
  }

  // h'kl = hkl - 0.5 sum_j (kj|jl)
  vector<double> buf(sizeij);
  int ij = 0;
  for (int i=0; i!=nocc; ++i) {
    for (int j=0; j<=i; ++j, ++ij) {
      buf[ij] = jop_->mo1e(j,i);
      for (int k=0; k!=nocc; ++k) {
        buf[ij] -= 0.5*first[(k+j*nocc)*mm+(k+i*nocc)];
      }
    }
  }
  copy(buf.begin(), buf.end(), jop_->mo1e().begin());
  delete[] first;

}

#endif
