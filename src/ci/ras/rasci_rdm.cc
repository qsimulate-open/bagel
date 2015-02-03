

#include <src/ci/ras/rasci.h>
//#include <src/util/prim_op.h>
//#include <src/util/math/algo.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;

void RASCI::compute_rdm12(shared_ptr<RASCivec> cbra, shared_ptr<RASCivec> cket) {
  cout << "compute_rdm12.." << endl;

  rdm1_ = make_shared<RDM<1>>(norb_);
  rdm2_ = make_shared<RDM<2>>(norb_);

  // since we consider here number conserving operators...
  auto dbra = make_shared<RASDvec>(cbra->det(), norb_*norb_);
  for (int ij = 0; ij != norb_*norb_; ++ij) {
    dbra->data(ij)->zero();
  }
//dbra->zero();
//sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<RASDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb_*norb_);
    for (int ij = 0; ij != norb_*norb_; ++ij) {
      dket->data(ij)->zero();
    }
  //dket->zero();
  //sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  tie(rdm1,rdm2) = compute_rdm12_last_step(dbra,dket,cbra);
  cout << "RASCI: RDM1" << endl;
  auto mat = rdm1->rdm1_mat(0);
  mat->print("1RDM",norb_);
//rdm1->print(1.0e-6);

//return compute_rdm12_last_step(dbra, dket, cbra);
}


void RASCI::sigma_2a1(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  cout << "sigma_2a1" << endl;
  const int nij = norb_*norb_;
  shared_ptr<const RASDeterminants> det = cc->det();

  for (auto& ispace : *det->stringspaceb()) {
    cout << "ispace" << endl;
    for (auto istring = ispace->begin(); istring != ispace->end(); ++istring) {
      cout << " istring" << endl;
      const std::bitset<nbit__> bbit = *istring;
      for (int ij = 0; ij != nij; ++ij) {    
        cout << "  ij" << endl;
        for (auto& phiblock : det->uncompressed_phia_ij(ij)) {
          cout << "   phiblock" << endl;
          for (auto& phi : phiblock) {
            cout << "    phi" << endl;
            const double sign = static_cast<double>(phi.sign);
            const auto sbit = det->string_bits_a(phi.source);
            cout << "     sbit:" << print_bit(sbit,norb_) << endl;
            const auto tbit = det->string_bits_a(phi.target);
            cout << "     tbit:" << print_bit(tbit,norb_) << endl;
          //if(!allowed(sbit,bbit)) continue; //alpha first
            if(!det->allowed(tbit,bbit)) continue;
            if(!det->allowed(sbit,bbit)) continue;
            cout << "     cc: " << cc->element(bbit,tbit) <<  endl;
            cout << "      d: " << d->data(ij)->element(bbit,sbit) << endl;
            d->data(ij)->element(bbit,sbit) += sign * cc->element(bbit,tbit); //beta first
            cout << "     ---: " << endl;
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
      const std::bitset<nbit__> abit = ispace->strings(ia);

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
      const std::bitset<nbit__> abit = *istring;
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
  RASCI::compute_rdm12_last_step(shared_ptr<RASDvec> dbra, shared_ptr<RASDvec> dket, shared_ptr<const RASCivec> cibra) const {

  cout << "last-step test" << endl;
  auto rdm2 = make_shared<RDM<2>>(norb_);
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


/*
  // 2RDM
  // \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm2 = make_shared<RDM<2>>(norb_);
  dgemm_("T", "N", ij, ij, nri, 1.0, dbra->data(0)->data(), nri, dket->data(0)->data(), nri, 0.0, rdm2->data(), ij);

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
*/
  return tie(rdm1, rdm2);
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
            if (std::abs(*i) > thr)
              tmp.emplace(-std::abs(*i), std::make_tuple(*i, ia, ib));
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

