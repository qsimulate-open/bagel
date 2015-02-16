

#include <src/ci/ras/rasci.h>
#include <src/util/prim_op.h>
//#include <src/util/math/algo.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;


void RASCI::compute_rdm12() {

  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);

  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_[0].
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  }
  if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
//auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, /*compressed=*/false, /*mute=*/true);
//cc_->set_det(detex);

//for (int i = 0; i != nstate_; ++i) compute_rdm12(i);

//cc_->set_det(det_);

  int istate = 0;
  tie(rdm1_av_,rdm2_av_) = compute_rdm12(cc_->data(istate), cc_->data(istate));
  rdm1_[0] = rdm1_av_;
  rdm2_[0] = rdm2_av_;
}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
     RASCI::compute_rdm12(shared_ptr<RASCivec> cbra, shared_ptr<RASCivec> cket) {
  cout << "compute_rdm12.." << endl;

//rdm1_ = make_shared<RDM<1>>(norb_);
//rdm2_ = make_shared<RDM<2>>(norb_);

  // since we consider here number conserving operators...
  auto dbra = make_shared<RASDvec>(cbra->det(), norb_*norb_);
  for (int ij = 0; ij != norb_*norb_; ++ij) {
    dbra->data(ij)->zero();
  }
  //new
  auto eket = make_shared<RASDvec>(cket->det(), norb_*norb_*norb_*norb_);
  for (int ij = 0; ij != norb_*norb_*norb_*norb_; ++ij)
    eket->data(ij)->zero();

//dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  sigma_2a1_new(cket, eket);
  sigma_2a2_new(cket, eket);
  sigma_2a3_new(cket, eket);
  sigma_2a4_new(cket, eket);

  shared_ptr<RASDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb_*norb_);
    for (int ij = 0; ij != norb_*norb_; ++ij) {
      dket->data(ij)->zero();
    }
  //dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  tie(rdm1,rdm2) = compute_rdm12_last_step(dbra,dket,cbra,eket);
  cout << "RASCI: RDM1" << endl;
  auto mat = rdm1->rdm1_mat(0);
  mat->print("1RDM",norb_);
//rdm1->print(1.0e-6);

  //RDM2 symmetrize (out-of-excitation free parts are copied)
  {
    for (int i = 0, ij = 0; i != norb_; ++i) {
      for (int j = 0; j != norb_; ++j, ++ij) {
        for (int k = 0, kl = 0; k != norb_; ++k) {
          for (int l = 0; l != norb_; ++l, ++kl) {
            if (kl > ij) {
              rdm2->element(i,j,k,l) = rdm2->element(k,l,i,j);
            }
          }
        }
      }
    }
  }

  //Trace
  { 
    double sum = 0.0;
    for (int j = 0; j != norb_; ++j) 
      sum += rdm1->element(j,j);
    cout << "1RDM Trace = " << sum << endl;
    sum = 0.0;
    for (int i = 0; i != norb_; ++i)
      for (int j = 0; j != norb_; ++j) 
        sum += rdm2->element(i,i,j,j);
    cout << "2RDM Trace = " << sum << endl;
  }
  //Partial trace
  {
    cout << "2RDM Partial Trace Sum_k (i,j,k,k)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1);
    for (int i = 0; i != norb_; ++i)
    for (int j = 0; j != norb_; ++j)
    for (int k = 0; k != norb_; ++k) {
      debug->element(i,j) -= 1.0/(nelea_+neleb_-1) * rdm2->element(i,j,k,k);
    }
    debug->print(1.0e-8);
  }
  {
    cout << "2RDM Partial Trace Sum_k (k,k,i,j)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1);
    for (int i = 0; i != norb_; ++i)
    for (int j = 0; j != norb_; ++j)
    for (int k = 0; k != norb_; ++k) {
      debug->element(i,j) -= 1.0/(nelea_+neleb_-1) * rdm2->element(k,k,i,j);
    }
    debug->print(1.0e-8);
  }

  //Symmetry check
  {//RDM1
    for (int j = 0; j != norb_; ++j)
    for (int i = 0; i != norb_; ++i) {
      double ij = rdm1->element(i,j);
      double ji = rdm1->element(j,i);
      if( abs(ij-ji) > 1.0e-10) { //assert(false); //cout << "ERROR1" << endl;
        cout << "R1: " << ij << " " << ji << endl;
      }
    }
  }
  {//RDM2 
    for (int l = 0; l != norb_; ++l)
    for (int k = 0; k != norb_; ++k)
    for (int j = 0; j != norb_; ++j)
    for (int i = 0; i != norb_; ++i) {
      double ijkl = rdm2->element(i,j,k,l);
      double klij = rdm2->element(k,l,i,j);
      double jilk = rdm2->element(j,i,l,k);
      double lkji = rdm2->element(l,k,j,i);
      if( abs(ijkl-klij) > 1.0e-10) { //assert(false); //cout << "ERROR1" << endl;
        cout << "W1: " << i << j << k << l << ":" << ijkl << " /= "  
                       << k << l << i << j << ":" << klij << endl;
      }
      if( abs(ijkl-jilk) > 1.0e-10) { // assert(false); //cout << "ERROR2" << endl;
        cout << "W2: " << i << j << k << l << ":" << ijkl << " /= "
                       << j << i << l << k << ":" << jilk << endl;
      }
      if( abs(ijkl-lkji) > 1.0e-10) { // assert(false); //cout << "ERROR3" << endl;
        cout << "W3: " << i << j << k << l << ":" << ijkl << " " << lkji << endl;
      }
    }
  }

  {//RDM2 print
    for (int l = 0; l != norb_; ++l)
    for (int k = 0; k != norb_; ++k)
    for (int j = 0; j != norb_; ++j)
    for (int i = 0; i != norb_; ++i) {
      double ijkl = rdm2->element(i,j,k,l);
      cout << "2RDM: " << i << j << k << l << ":" << ijkl << endl;
    }
  }

  //Energy calculation
  cout << "RASCI: Energy calculated from RDM:" << endl;
  cout << "Number of closed orbitals: " << ncore_ << endl;
  cout << "Number of active orbitals: " << norb_  << endl;
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
  cout << "Nuc       = " << geom_->nuclear_repulsion() << endl;
  cout << "Core      = " << jop_->core_energy() << endl;
  cout << "Nuc+Core  = " << nuc_core << endl;

  shared_ptr<const Matrix> h1 = jop_->mo1e()->matrix();
  assert(norb_ == h1->ndim() && norb_ == h1->mdim());
  h1->print("1e integral",norb_);
  auto onerdm = rdm1->rdm1_mat(0);
  double  e1 = ddot_(norb_*norb_, h1->element_ptr(0,0), 1, onerdm->element_ptr(0,0), 1);
  cout << "1E energy = " << e1 << endl;

  shared_ptr<const Matrix> pint2 = jop_->mo2e()->matrix();
  auto int2 = make_shared<Matrix>(norb_*norb_*norb_*norb_,1);
  sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), norb_, norb_, norb_, norb_); //conver to chemist not.

  auto low = {0,0,0,0};
  auto up  = {norb_,norb_,norb_,norb_};
  auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); 
  auto twordm = make_shared<Matrix>(norb_*norb_*norb_,norb_,1); 
  copy(view.begin(), view.end(), twordm->begin());

  double e2 = 0.5 * ddot_(norb_*norb_*norb_*norb_, int2->element_ptr(0,0), 1, twordm->element_ptr(0,0), 1);
  cout << "2E energy = " << e2 << endl;

  //Energy print
  cout << "Total energy = " << nuc_core + e1 + e2 << endl;

//return compute_rdm12_last_step(dbra, dket, cbra);
  return tie(rdm1, rdm2);
}


void RASCI::sigma_2a1(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a1" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) { 
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
    //const bitset<nbit__> bbit = ispace->strings(ib);
      const auto bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) { 
          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_a(phi.source);
            const auto tbit = det->string_bits_a(phi.target);
            const auto ij = phi.ij;
            if(!det->allowed(tbit,bbit)) continue;
            if(!det->allowed(sbit,bbit)) continue;
            d->data(ij)->element(bbit,sbit) += sign * cc->element(bbit,tbit);
          }
        }
      }
    }
  }

}

void RASCI::sigma_2a2(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a2" << endl;
  
//const int nij = norb_*norb_;
  shared_ptr<const RASDeterminants> det = cc->det();
//const size_t lb = det->lenb();

  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
  //const size_t offset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) { // determinants associated with a given space
    //const bitset<nbit__> abit = ispace->strings(ia);
      const auto abit = ispace->strings(ia);

    //for (auto& phi : det->uncompressed_phib(ia+offset)) {
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { // determinants associated with a given space

          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source);
            const auto tbit = det->string_bits_b(phi.target);
            const auto ij = phi.ij;
          //if(i == j) {
          //  cout << "diag[" << i << "] : " << phi.source << " " << phi.target << endl;
          //}
            if(!det->allowed(abit,tbit)) continue;
            if(!det->allowed(abit,sbit)) continue;
            d->data(ij)->element(sbit,abit) += sign * cc->element(tbit,abit);
          }
        }
      }
    }
  }
}
/*
  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
    for (auto istring = ispace->begin(); istring != ispace->end(); ++istring) {
      const bitset<nbit__> abit = *istring;
    //for (int ij = 0; ij != nij; ++ij) {    
      for (int i = 0, ij = 0; i != norb_; ++i)
      for (int j = 0; j !=norb_; ++j, ++ij) {
        for (auto& phiblock : det->uncompressed_phib_ij(ij)) { //<DetMapBlock>
          for (auto& phi : phiblock) { //<DetMap>
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source);
            const auto tbit = det->string_bits_b(phi.target);
            if(i == j) {
              cout << "diag[" << i << "] : " << phi.source << " " << phi.target << endl;
            }
            if(!det->allowed(abit,tbit)) continue;
            if(!det->allowed(abit,sbit)) continue;
            d->data(ij)->element(sbit,abit) += sign * cc->element(tbit,abit);
          }
        }
      }
    }
  }
}
*/
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  RASCI::compute_rdm12_last_step(shared_ptr<RASDvec> dbra, shared_ptr<RASDvec> dket, shared_ptr<const RASCivec> cibra, shared_ptr<RASDvec> eket) const {

  cout << "last-step entered.." << endl;
/*
  const int nri = dbra->lena()*dbra->lenb();
  const int ij  = norb_*norb_;

  if (nri != dket->lena()*dket->lenb())
    throw logic_error("FCI::compute_rdm12_last_step called with inconsistent RI spaces");
*/

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  shared_ptr<const RASDeterminants> det = cibra->det();
  auto rdm1 = make_shared<RDM<1>>(norb_);
  rdm1->zero();
//dgemv_("T", nri, ij, 1.0, dket->data(0)->data(), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);
  for (int i = 0, ij = 0; i != norb_; ++i){
    for (int j = 0; j!= norb_; ++j, ++ij) {

      for (auto& cblock : cibra->blocks()) {
        if (!cblock) continue;
  
        for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
          for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
          //auto cptr = cblock->data() + cab;
          //cout << "[" << ca << "," << cb << "] = " << *cptr << " " <<cblock->element(cab) << endl;
            auto abit = cblock->stringsa()->strings(ca);
            auto bbit = cblock->stringsb()->strings(cb);
          //rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cblock->element(cab);
          //rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cblock->element(bbit,abit); 
            rdm1->element(i,j) += dket->data(ij)->element(bbit,abit) * cibra->element(bbit,abit); 
          //cout << endl;
          //cout << "check: " << cblock->element(cab) << endl;
          //cout << "       " << cblock->element(bbit,abit) << endl;
          //cout << "       " << cibra->element(bbit,abit) << endl;
          //cout << endl;
          }
        }
      }
    }
  }
  cout << "1rdm done.." << endl;

  auto new2 = make_shared<RDM<2>>(norb_);
  {
    //NEW 2RDM
    new2->zero();
    for (int i = 0, ij = 0, ijkl = 0; i != norb_; ++i){
      for (int j = 0; j != norb_; ++j, ++ij) {
        //fixed ij
        for (int k = 0, kl = 0; k != norb_; ++k){
          for (int l = 0; l != norb_; ++l, ++kl, ++ijkl) {
            //fixed kl
            
            for (auto& cblock : cibra->blocks()) { //just run over allowed blocks/ TODO: replace this
              if (!cblock) continue;
      
              for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
                for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
                  auto abit = cblock->stringsa()->strings(ca);
                  auto bbit = cblock->stringsb()->strings(cb);
                  assert(det->allowed(abit,bbit));
                  new2->element(k,l,i,j) += cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit);
                //new2->element(i,j,k,l) += cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit);
                //if(abs(cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit)) > 0.01) {
                //}
                }
              }
 
            }
          }
        }
 
      }
    }
    // put in diagonal into 2RDM
    // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
    for (int i = 0; i != norb_; ++i) {
      for (int k = 0; k != norb_; ++k) {
        for (int j = 0; j != norb_; ++j) {
          new2->element(j,k,k,i) -= rdm1->element(j,i);
        } 
      }
    }

  cout << "2rdm done.." << endl;
 
  }
/*
  // 2RDM
  // \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm2 = make_shared<RDM<2>>(norb_);
  rdm2->zero();
//dgemm_("T", "N", ij, ij, nri, 1.0, dbra->data(0)->data(), nri, dket->data(0)->data(), nri, 0.0, rdm2->data(), ij);
  for (int i = 0, ij = 0; i != norb_; ++i){
    for (int j = 0; j != norb_; ++j, ++ij) {
      //fixed ij
      for (int k = 0, kl = 0; k != norb_; ++k){
        for (int l = 0; l != norb_; ++l, ++kl) {
          //fixed kl
          
          for (auto& cblock : cibra->blocks()) { //just run over allowed blocks/ TODO: replace this
            if (!cblock) continue;
    
            for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
              for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
                auto abit = cblock->stringsa()->strings(ca);
                auto bbit = cblock->stringsb()->strings(cb);
                rdm2->element(i,j,k,l) += 2.0* dbra->data(ij)->element(bbit,abit) * dket->data(kl)->element(bbit,abit); //ab & ba
              }
            }

          }
        }
      }

    }
  }

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, &rdm2->element(0,0,k,i));
    }
  }
  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb_; ++i)
    for (int k = 0; k != norb_; ++k)
      for (int j = 0; j != norb_; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  //add
  daxpy_(norb_*norb_*norb_*norb_, 1.0, &rdm2->element(0,0,0,0), 1, &new2->element(0,0,0,0), 1);

*/

//return tie(rdm1, rdm2);
  return tie(rdm1, new2);
}

void RASCI::sigma_2a1_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a1_new" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) { 
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) { 
          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_a(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(tbit,bbit)) continue;
            for (auto& phi2 : det->phia(phi.source)) {
              assert(phi2.target == phi.source);
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbit = det->string_bits_a(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sbit,bbit)) continue;
              d->data(ij + kl*norb_*norb_)->element(bbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}

void RASCI::sigma_2a2_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a2_new" << endl;
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspacea()) { 
    for (size_t ia = 0; ia != ispace->size(); ++ia) {
    //const bitset<nbit__> abit = ispace->strings(ia); //fix alpha-det
      const auto abit = ispace->strings(ia);
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset(); //scan through beta-det
        for (size_t jb = 0; jb != jspace->size(); ++jb) { 
          for (auto& phi : det->phib(jb+offset)) { //first E
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_b(phi.target); //coeff
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue; //coeff
            for (auto& phi2 : det->phib(phi.source)) {
              assert(phi2.target == phi.source);
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbit = det->string_bits_b(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(abit,sbit)) continue;
              d->data(ij + kl*norb_*norb_)->element(sbit,abit) += sign * sign2 * cc->element(tbit,abit);
            //if(ij == 5 && kl == 30) cout << "0550 : " << sign * sign2 * cc->element(tbit,abit) << endl;
            //if(ij == 30 && kl == 5) cout << "5005 : " << sign * sign2 * cc->element(tbit,abit) << endl;
            }
          }
        }
      }
    }
  }

  cout << "break" << endl;

}

void RASCI::sigma_2a3_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a3_new" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) { 
    const size_t boffset = ispace->offset();
    for (size_t ib = 0; ib != ispace->size(); ++ib) {
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) {  //fix a-string
          for (auto& phi : det->phia(ja+offset)) {
            assert(phi.target == ja+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_a(phi.source); //sbit(alpha)
            const auto tbit = det->string_bits_a(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(tbit,bbit)) continue; //coeff
            for (auto& phi2 : det->phib(ib+boffset)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbbit = det->string_bits_b(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sbit,sbbit)) continue; //double replacement ket
            //cout << cc->element(bbit,tbit) << endl;
            //cout << d->data(ij + kl*norb_*norb_)->element(sbbit,sbit) << endl;
            //cout << endl;
              d->data(ij + kl*norb_*norb_)->element(sbbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            //d->data(kl + ij*norb_*norb_)->element(sbbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}

void RASCI::sigma_2a4_new(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a4_new" << endl;
  //based on sigma_2a2
  
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspacea()) { 
    const size_t aoffset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) {
      const bitset<nbit__> abit = ispace->strings(ia);
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) {  //fix b-string
          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_b(phi.source); //sbit(beta)
            const auto tbit = det->string_bits_b(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue; //coeff
            for (auto& phi2 : det->phia(ia+aoffset)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sabit = det->string_bits_a(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(sabit,sbit)) continue; //double replacement ket
            //cout << cc->element(bbit,tbit) << endl;
            //cout << d->data(ij + kl*norb_*norb_)->element(sbbit,sbit) << endl;
            //cout << endl;
              d->data(ij + kl*norb_*norb_)->element(sbit,sabit) += sign * sign2 * cc->element(tbit,abit);
            //d->data(kl + ij*norb_*norb_)->element(sbbit,sbit) += sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

}


// note that this does not transform internal integrals (since it is not needed in CASSCF).
pair<shared_ptr<Matrix>, vector<double>> RASCI::natorb_convert() {
  assert(rdm1_av_ != nullptr);
  pair<shared_ptr<Matrix>, vector<double>> natorb = generate_natural_orbitals();
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  return natorb;
}


void RASCI::update_rdms(const shared_ptr<Matrix>& coeff) {
  (rdm1_av_->rdm1_mat(0))->print("ORIGINAL RDM1");
  for (auto iter = rdm1_.begin(); iter != rdm1_.end(); ++iter)
    (*iter)->transform(coeff);
  for (auto iter = rdm2_.begin(); iter != rdm2_.end(); ++iter)
    (*iter)->transform(coeff);
  (rdm1_av_->rdm1_mat(0))->print("TRANSFORMD RDM1");
/* TODO: 1state only at the moment
  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_.size() > 1)   rdm1_av_->transform(coeff);
  if (rdm2_.size() > 1)   rdm2_av_->transform(coeff);
  assert(rdm1_.size() > 1 || rdm1_.front() == rdm1_av_);
  assert(rdm2_.size() > 1 || rdm2_.front() == rdm2_av_);
*/


  //Energy calculation
  cout << "T// RASCI: Energy calculated from RDM:" << endl;
  cout << "T// Number of closed orbitals: " << ncore_ << endl;
  cout << "T// Number of active orbitals: " << norb_  << endl;
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
  cout << "T// Nuc       = " << geom_->nuclear_repulsion() << endl;
  cout << "T// Core      = " << jop_->core_energy() << endl;
  cout << "T// Nuc+Core  = " << nuc_core << endl;

  shared_ptr<const Matrix> h1 = jop_->mo1e()->matrix();
  assert(norb_ == h1->ndim() && norb_ == h1->mdim());
//h1->print("1e integral",norb_);
  auto onerdm = rdm1_av_->rdm1_mat(0);
  double  e1 = ddot_(norb_*norb_, h1->element_ptr(0,0), 1, onerdm->element_ptr(0,0), 1);
  cout << "T// 1E energy = " << e1 << endl;

  shared_ptr<const Matrix> pint2 = jop_->mo2e()->matrix();
  auto int2 = make_shared<Matrix>(norb_*norb_*norb_*norb_,1);
  sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), norb_, norb_, norb_, norb_); //conver to chemist not.

  auto low = {0,0,0,0};
  auto up  = {norb_,norb_,norb_,norb_};
  auto view = btas::make_view(rdm2_av_->range().slice(low,up), rdm2_av_->storage()); 
  auto twordm = make_shared<Matrix>(norb_*norb_*norb_,norb_,1); 
  copy(view.begin(), view.end(), twordm->begin());

  double e2 = 0.5 * ddot_(norb_*norb_*norb_*norb_, int2->element_ptr(0,0), 1, twordm->element_ptr(0,0), 1);
  cout << "T// 2E energy = " << e2 << endl;

  //Energy print
  cout << "T// Total energy = " << nuc_core + e1 + e2 << endl;

}

pair<shared_ptr<Matrix>, vector<double>> RASCI::generate_natural_orbitals() {
  cout << "GENERATE RAS NATURAL ORBITALS .... " << endl;
  shared_ptr<Matrix> rdm1 = rdm1_av_->rdm1_mat(/*closed*/0);
  //subspace rdm
  vector<shared_ptr<Matrix>> srdm;
  srdm.push_back(make_shared<Matrix>(*rdm1->get_submatrix(0,0,ras_[0],ras_[0])));
  srdm.push_back(make_shared<Matrix>(*rdm1->get_submatrix(ras_[0],ras_[0],ras_[1],ras_[1])));
  srdm.push_back(make_shared<Matrix>(*rdm1->get_submatrix(ras_[0]+ras_[1],ras_[0]+ras_[1],ras_[2],ras_[2])));

  auto outm = make_shared<Matrix>(norb_,norb_,true);
  outm->unit();
  vector<double> outv(norb_);

  for(int ispace = 0; ispace != 3; ++ispace) {
//for(int ispace = 1; ispace != 2; ++ispace) {

    auto buf = make_shared<Matrix>(ras_[ispace],ras_[ispace],true);
    buf->add_diag(2.0);
    daxpy_(ras_[ispace]*ras_[ispace], -1.0, srdm[ispace]->data(), 1, buf->data(), 1);
 
    VectorB vec(ras_[ispace]);
    buf->diagonalize(vec);
    for (int i = 0; i != ras_[ispace]; ++i) {
      cout << "vec (" << i << ") = " << vec(i) << endl;
    }
 
    for (auto& i : vec) i = 2.0-i;
 
    map<int,int> emap;
    auto buf2 = buf->clone();
    VectorB vec2(ras_[ispace]);
    // sort eigenvectors so that buf is close to a unit matrix
    // target column
    for (int i = 0; i != ras_[ispace]; ++i) {
      // first find the source column
      tuple<int, double> max = make_tuple(-1, 0.0);
      for (int j = 0; j != ras_[ispace]; ++j)
        if (fabs(buf->element(i,j)) > get<1>(max))
          max = make_tuple(j, fabs(buf->element(i,j)));
 
      // register to emap
      if (emap.find(get<0>(max)) != emap.end()) throw logic_error("this should not happen. RDM<1>::generate_natural_orbitals()");
      emap.emplace(get<0>(max), i);
 
      // copy to the target
      copy_n(buf->element_ptr(0,get<0>(max)), ras_[ispace], buf2->element_ptr(0,i));
      vec2(i) = vec(get<0>(max));
    }
 
    // fix the phase
    for (int i = 0; i != ras_[ispace]; ++i) {
      if (buf2->element(i,i) < 0.0)
        blas::scale_n(-1.0, buf2->element_ptr(0,i), ras_[ispace]);
    }
 
    for (int i = 0; i != ras_[ispace]; ++i) {
      cout << "vec2(" << i << ") = " << setw(10) << setprecision(6) << vec2(i) << endl;
    }

    assert(buf2->ndim() == ras_[ispace]);
    assert(buf2->mdim() == ras_[ispace]);
    assert(vec2.size() == ras_[ispace]);
    switch (ispace) {
      case 0: //RAS1
        outm->copy_block(0, 0, buf2->ndim(), buf2->mdim(), buf2);
        for (int i = 0; i != vec2.size(); ++i) {
          outv[i] = vec2(i);
        }
        break;
      case 1: //RAS2
        outm->copy_block(ras_[0], ras_[0], buf2->ndim(), buf2->mdim(), buf2);
        for (int i = 0; i != vec2.size(); ++i) {
          outv[ras_[0]+i] = vec2(i);
        }
        break;
      case 2: //RAS3
        outm->copy_block(ras_[0]+ras_[1], ras_[0]+ras_[1], buf2->ndim(), buf2->mdim(), buf2);
        for (int i = 0; i != vec2.size(); ++i) {
          outv[ras_[0]+ras_[1]+i] = vec2(i);
        }
        break;
      default:
        cout << "out of RAS space" << endl;
        assert(false);
        break;
    }
  } //ispace

  for (int i = 0; i != norb_; ++i) {
    cout << "nat occ(" << i << ") = " << setw(10) << setprecision(6) << outv[i] << endl;
  }
  
  outm->purify_unitary();
  outm->print("TRANSFORMATION MATRIX");

  return {outm, outv}; //TEMP
}


#if 0

  shared_ptr<const RASDeterminants> det = cc->det();
  const int norb = det->norb();
  const size_t lb = det->lenb();

  const double* const source_base = cc->data();

  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) { //E_{i'j}

      double* const target_base = d->data(ij)->data();

      for (auto& target_bspace : *det->stringspaceb()) {
        const size_t tlb = target_bspace->size();
     
        // looping over source_aspace
        for (auto& phiblock : det->phia_ij(ij) ) {
          const shared_ptr<const RASString>& source_aspace = phiblock.source_space();
     
          for (auto& phi : phiblock) {
          //auto target_aspace = det->space<0>(det->string_bits_a(phi.target));
            const double sign = static_cast<double>(phi.sign);
            double* const target_array = target_base + phi.source*tlb;
            cout << "SIGN: " << sign << endl;
            cout << "print " << tlb << lb << &source_base << &target_array << endl;
          }
        }
      }
    }
  }
    



  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
//  const size_t offset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) { // determinants associated with a given space
//    const double* const source_array0 = cc->element_ptr(0,offset+ia); //pointer to coeff array associated with determinant |ia>
    
      for (int ij = 0; ij != nij; ++ij) { //E_ij (i>=j or i<=j)
//      double* const target_array0 = d->data(ij)->element_ptr(0,offset+ia);
        for (auto& phiblock : det->phib_ij(ij)) { //allowed beta excitation spaces
          for (auto& phi : phiblock) { //single determinant 
            const double sign = static_cast<double>(phi.sign); 
          //d->data(ij)->element(phi,ia) += sign   
          //target_array0[phi.source] += sign * source_array0[phi.target];
          //target_array0[phi.source] += sign * source_array0[phi.target];
          }
        }
      }
    }
  }

  //coeff blocks (target)
  for (auto& cblock : cc->blocks()) {
    if (!cblock) continue;
    auto c_aspace = cblock->stringsa(); //a-space of coff. block
    auto c_bspace = cblock->stringsb(); //b-space of  "

    for (int ij = 0; ij !=nij; ++ij) { //Dvec(ij)
      for (auto& phiblock : det->target_phib_ij(ij) ) { 
        if(!c_bspace->matches(phiblock.target_space())) continue; // coeff beta space must match with target space of displacement
    
        for (auto& phi : phiblock) { //<DetMap>
          auto sdet = det->string_bits_b(phi.source);
          auto tdet = det->string_bits_b(phi.target);
        }


        //<DetMapBlock>
        const shared_ptr<const RASString>& source_bspace = phiblock.source_space();
        vector<tuple</*source*/size_t,/*sign*/int,/*offset_to_target*/size_t> reduced_phi;
        for (auto& phi : phiblock) { 
          //<DetMap>
          auto target_bspace = det->space<1>(det->string_bits_b(phi.target));
          if (det->allowed(source_aspace, target_bstace)){
            
          }
        }

      }// phiblock      


    } //cblock    

  }//ij


  size_t offset = 0;
  for (auto& ispace : *det->stringspacea()) {
    //fill C for given string space
    const size_t la = ispace->size(); //#of a-strings
    Matrix C(lb,la);
    C.zero();
    cout << "Column matrix [" << C.ndim() << "," << C.mdim() << "]" << endl;
  //double* const cdata = C.element_ptr(0,ia);
    for (auto& cblock : cc->blocks()) {
      if (!cblock) continue;
      if (!ispace->matches(cblock->stringsa())) continue;

      //matching a-block
      //check b-block and copy to C
      for (size_t ca = 0, cab = 0; ca < cblock->stringsa()->size(); ++ca) {
        for (size_t cb = 0; cb < cblock->stringsb()->size(); ++cb, ++cab) {
          auto cptr = cblock->data() + cab;
          cout << "[" << ca << "," << cb << "] = " << *cptr << " " <<cblock->element(cab) << endl;
          cout << "A" << print_bit(cblock->stringsa()->strings(ca),norb_) << " : " << endl;
          cout << "B" << print_bit(cblock->stringsb()->strings(cb),norb_) << " : " << endl;
          //resolve ca, cb and get lexical number
    


        }
      }
      offset += la;
    }
    
    //fill D
//Matrix D(lb,nij);
//  D.zero();

//  for (size_t ia = 0; ia < ispace->size(); ++ia) {
//    cout << ia << endl;
//  }//ia

  }//ispace

  for (int ij = 0; ij != nij; ++ij) {
    cout << "IJ = " << ij << endl;
    for (auto& jblock : cc->det()->uncompressed_phib_ij(ij)) {
    //b <DetMapBlock>
      cout << "  space =" << jblock.source_space()->nholes() << " " << jblock.source_space()->nparticles() << endl;
      for (auto& jb : jblock) { //DetMap
        DetMap temp = jb;
        cout << "       " <<  temp.ij << endl; // this is same as ij above
      }
    }
  }
  assert(false);

}

  const size_t lb = det->lenb();
  VectorB T(lb); //target: Dvec
  VectorB S(lb); //source: coeff

  for (auto& ispace : *det->stringspacea()) { //a-space
    for (size_t ia = 0; ia < ispace->size()) { //a-string
      S.zero();
      //fill S
      for (auto& iblock : cc->block()) {
        if (!iblock) continue;
      }
      T.zero();
  //Source: coeff
  for (auto& iblock : cc->blocks()) { // blocks of RAS coeff
    if (!iblock) continue; // allowed blocks by RAS def 
//  double* source_base = iblock->data(); // pointer to the block
    for (auto& ia : *iblock->stringsa()) { //a-strings in iblock (fixed)
    //single a-string <DetMap>
      for (int ij = 0; ij != nij; ++ij) { //component of RASDvec
//      double* target_base = d->data(ij);

        //target column in Dvec[ij]
        T.zero()
        for (auto& jblock : cc->det()->uncompressed_phib_ij(ij)) { //b-string vectors with fixed E_ij
        //beta <DetMapBlock>
          if(!det->allowed(iblock->stringsa(),jblock->source_space())) continue;
          //            RASString
          for(auto& jb : jblock) { //b-string
            //<DetMap>
            const double sign = static_cast<double>(jb.sign);
            T[jb.source] += sign * source [jb.target]
          }

        }
      }
    }
  }

  //const DataType* element_ptr(size_t i, size_t j) const { return cc()+i+j*lenb_; }
    int icounter = 0
    for (auto& ia : *iblock->stringsa()) { //alpha determinants associated with the given block
      const double* const source_array0 = source_base + ?*lenb;



        double* i = iblock->data();
        for (auto& ia : *iblock->stringsa()) {
          for (auto& ib : *iblock->stringsb()) {
            if (abs(*i) > thr)
              tmp.emplace(-abs(*i), make_tuple(*i, ia, ib));
            ++i;
          }
        }
      }
      for (auto& i : tmp)


      for (int ij = 0; ij != nij; ++ij) {
        double* const target_array0 = d->data(ij)

        for (auto& iter : cc->det()->phib(ij)) {
          const double sign = static_cast<double>(iter.sign);
          target_array0[
        }
      }
    }

  for (auto& ispace : *det->stringspacea()) { // alpha determinant space
    const size_t offset = ispace->offset();
    for (size_t ia = 0; ia != ispace->size(); ++ia) { // determinants associated with a given space
//    const double* const source_array0 = cc->element_ptr(0,offset+ia); //pointer to coeff array associated with determinant |ia>
    
      for (int ij = 0; ij != nij; ++ij) { //E_ij (i>=j or i<=j)
//      double* const target_array0 = d->data(ij)->element_ptr(0,offset+ia);
        for (auto& phiblock : det->phib_ij(ij)) { //allowed beta excitation spaces
          for (auto& phi : phiblock) { //single determinant 
            const double sign = static_cast<double>(phi.sign); 
          //d->data(ij)->element(phi,ia) += sign   
          //target_array0[phi.source] += sign * source_array0[phi.target];
          //target_array0[phi.source] += sign * source_array0[phi.target];
          }
        }
      }
    }
  }

#endif

