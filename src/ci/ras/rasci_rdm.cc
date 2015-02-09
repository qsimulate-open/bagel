

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
      const bitset<nbit__> bbit = ispace->strings(ib);
      for (auto& jspace : *det->stringspacea()) {
        const size_t offset = jspace->offset();
        for (size_t ja = 0; ja != jspace->size(); ++ja) { 
          for (auto& phi : det->phia(ja+offset)) {
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
      const bitset<nbit__> abit = ispace->strings(ia);

    //for (auto& phi : det->uncompressed_phib(ia+offset)) {
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { // determinants associated with a given space

          for (auto& phi : det->phib(jb+offset)) {
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
                  new2->element(i,j,k,l) += cibra->element(bbit,abit) * eket->data(ijkl)->element(bbit,abit);
                }
              }
 
            }
          }
        }
 
      }
    }
 
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
    // put in diagonal into 2RDM
    // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
    for (int i = 0; i != norb_; ++i)
      for (int k = 0; k != norb_; ++k)
        for (int j = 0; j != norb_; ++j)
          new2->element(j,k,k,i) -= rdm1->element(j,i);

  cout << "2rdm done.." << endl;

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
      const bitset<nbit__> abit = ispace->strings(ia);
      for (auto& jspace : *det->stringspaceb()) {
        const size_t offset = jspace->offset();
        for (size_t jb = 0; jb != jspace->size(); ++jb) { 
          for (auto& phi : det->phib(jb+offset)) {
            assert(phi.target == jb+offset);
            const double sign = static_cast<double>(phi.sign);
          //const auto ibit = det->string_bits_b(phi.source); //intermediate
            const auto tbit = det->string_bits_b(phi.target);
            const auto kl = phi.ij;
            if(!det->allowed(abit,tbit)) continue;
            for (auto& phi2 : det->phib(phi.source)) {
              const double sign2 = static_cast<double>(phi2.sign);
              const auto sbit = det->string_bits_b(phi2.source);
              const auto ij = phi2.ij;
              if(!det->allowed(abit,sbit)) continue;
              d->data(ij + kl*norb_*norb_)->element(sbit,abit) += sign * sign2 * cc->element(tbit,abit);
            }
          }
        }
      }
    }
  }

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
              d->data(ij + kl*norb_*norb_)->element(sbbit,sbit) += 2.0 * sign * sign2 * cc->element(bbit,tbit);
            }
          }
        }
      }
    }
  }

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

