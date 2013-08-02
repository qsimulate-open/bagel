//
// BAGEL - Parallel electron correlation program.
// Filename: relmofile.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


//plan is to replace zmofile with relmofile because zmofile does not have complex integrals anyway

RelMOFile::RelMOFile(const shared_ptr<const Reference> ref, const string method)
: geom_(ref->geom()), ref_(ref), core_fock_(make_shared<Matrix>(geom_->nbasis(), geom_->nbasis())), coeff_(ref_->coeff()) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}

RelMOFile::RelMOFile(const shared_ptr<const Reference> ref, const shared_ptr<const Coeff> c, const string method)
: hz_(false), geom_(ref->geom()), ref_(ref), core_fock_(make_shared<Matrix>(geom_->nbasis(), geom_->nbasis())), coeff_(c) {

  do_df_ = geom_->df().get();
  if (!do_df_) throw runtime_error("for the time being I gave up maintaining non-DF codes.");

  hz_ = (method=="HZ");
}

RelMOFile::~RelMOFile() {
}

//one line change in order to add a function to block diagonalize
double ZMOFile::create_Jiiii(const int nstart, const int nfence) {
  // first compute all the AO integrals in core
  nocc_ = nfence - nstart;
  nbasis_ = geom_->nbasis();
  const int nbasis = nbasis_;

  // one electron part
  double core_energy;
  shared_ptr<const ZMatrix> buf1e;
  tie(buf1e, core_energy) = compute_mo1e(nstart, nfence);

  // two electron part.
  // this fills mo2e_1ext_ and returns buf2e which is an ii/ii quantity
  unique_ptr<complex<double>[]> buf2e = compute_mo2e(nstart, nfence);

  //block diaganolize according kramers symmetry
  //
  //

  compress(buf1e, buf2e);
  return core_energy;
}

//
//Function for block diagonalization of integrals
//

//does not need to be changed as long as product of block diagonalization is same type and format
void ZMOFile::compress(shared_ptr<const ZMatrix> buf1e, unique_ptr<complex<double>[]>& buf2e) {

  const int nocc = nocc_;
  sizeij_ = nocc*nocc;
  mo2e_ = unique_ptr<complex<double>[]>(new complex<double>[sizeij_*sizeij_]);
  copy_n(buf2e.get(), sizeij_*sizeij_, mo2e_.get());

  // h'ij = hij - 0.5 sum_k (ik|kj)
  const int size1e = nocc*nocc;
  mo1e_ = unique_ptr<complex<double>[]>(new complex<double>[size1e]);
  int ij = 0;
  for (int i = 0; i != nocc; ++i) {
    for (int j = 0; j != nocc; ++j, ++ij) {
      mo1e_[ij] = buf1e->element(j,i);
      for (int k = 0; k != nocc; ++k)
        mo1e_[ij] -= 0.5*buf2e[j+k*nocc+k*nocc*nocc+i*nocc*nocc*nocc];
    }
  }
}

//
//Relativistic 1 and 2e functions in progress
//
//inputs?
//ZJOP or ZMOFILE needs Geometry geom and Reference ref
tuple<shared_ptr<const ZMatrix>, double> ZJop::compute_mo1e(const int nstart, const int nfence) {
  const size_t nbasis = geom_->nbasis();

  double core_energy = 0;

  auto geom_ = geom->relativistic( "FALSE" );
  auto relhcore = make_shared<RelHCore>(geom_);
  auto relcoeff = ref->relcoeff();

  auto tmp = unique_ptr<complex<double>[]>(new complex<double>[8*nbasis*nbasis]);
  auto out = make_shared<ZMatrix>(2*nbasis, 2*nbasis);
  //TODO conjugated?
  // Hij = relcoeff(T) * relhcore * relcoeff
  zgemm3m(c, n, 2*nbasis, 4*nbasis, 4*nbasis, relcoeff->data(), 2*nbasis, relhcore->data(), 4*nbasis, tmp->get(), 2*nbasis);
  zgemm3m(n, n, 2*nbasis, 4*nbasis, 4*nbasis, tmp->get(), 2*nbasis, relcoeff->data(), 4*nbasis, out->data(), 2*nbasis);

  //TODO include some density adjustment? see zmofile
  core_energy = relhcore.trace();

  return make_tuple(out, core_energy);
}

//inputs?
unique_ptr<complex<double>[]> ZJop::compute_mo2e(const int nstart, const int nfence) {

//slightly modified code from rel/dmp2.cc to form 3 index integrals that we can build into 4 index with form4index
  const size_t nbasis = geom_->nbasis();
  shared_ptr<const RelReference> ref = dynamic_pointer_cast<const RelReference>(ref_);

  assert(nbasis*4 == ref->relcoeff()->ndim());
  assert(nbasis*2 == ref->relcoeff()->mdim());

  // Separate Coefficients into real and imaginary
  // correlated occupied orbitals
  array<shared_ptr<const Matrix>, 4> rocoeff;
  array<shared_ptr<const Matrix>, 4> iocoeff;
  for (int i = 0; i != 4; ++i) {
    shared_ptr<const ZMatrix> oc = ref->relcoeff()->get_submatrix(i*nbasis, 0, nbasis, 4*nbasis);
    rocoeff[i] = oc->get_real_part();
    iocoeff[i] = oc->get_imag_part();
  }

  // (1) make DFDists
  shared_ptr<Geometry> cgeom;
  vector<shared_ptr<const DFDist>> dfs;
  if (abasis_.empty()) {
    dfs = geom_->dfs()->split_blocks();
    dfs.push_back(geom_->df());
  } else {
    auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
    cgeom = make_shared<Geometry>(*geom_, info, false);
    cgeom->relativistic(false /* do_gaunt */);
    dfs = cgeom->dfs()->split_blocks();
    dfs.push_back(cgeom->df());
  }
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  // (2) first-transform
  list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, rocoeff, iocoeff);

  // (3) split and factorize
  list<shared_ptr<RelDFHalf>> half_complex_exch;
  for (auto& i : half_complex) {
    list<shared_ptr<RelDFHalf>> tmp = i->split(false);
    half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
  }
  half_complex.clear();
  DFock::factorize(half_complex_exch);

  // (4) compute (gamma|ia)
  list<shared_ptr<RelDFFull>> dffull;
  for (auto& i : half_complex_exch)
    dffull.push_back(make_shared<RelDFFull>(i, rvcoeff, ivcoeff));
  DFock::factorize(dffull);
  dffull.front()->scale(dffull.front()->fac()); // take care of the factor
  assert(dffull.size() == 1);
  shared_ptr<const RelDFFull> full = dffull.front();
//end mp2 copied code

  // use form_4index function to product 4 index (ij|kl) = sum (ij|gamma)(gamma|kl)


}

