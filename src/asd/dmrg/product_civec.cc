//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: product_civec.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/asd/dmrg/product_civec.h>
#include <src/ci/ras/civec_spinop.h>

using namespace std;
using namespace bagel;

RASBlockVectors RASBlockVectors::transpose_civecs(shared_ptr<const RASDeterminants> transdet) const {
  if (!transdet) transdet = det()->transpose();
  const int M = mdim();
  const int nasign = 1 - (((det()->nelea()*det()->neleb())%2) << 1);

  RASBlockVectors out(transdet, M);

  for (auto& source_block : det()->blockinfo()) {
    if (!source_block->empty()) {
      auto target_block = transdet->blockinfo(source_block->stringsa(), source_block->stringsb());
      assert(!target_block->empty());
      for (int m = 0; m < M; ++m) {
        const double* sourcedata = element_ptr(source_block->offset(), m);
        double* targetdata = out.element_ptr(target_block->offset(), m);
        blas::transpose(sourcedata, source_block->lenb(), source_block->lena(), targetdata, static_cast<double>(nasign));
      }
    }
  }

  return out;
}


// Constructor
ProductRASCivec::ProductRASCivec(shared_ptr<RASSpace> space, shared_ptr<const DMRG_Block> left, const int nelea, const int neleb) :
  space_(space), left_(left), nelea_(nelea), neleb_(neleb) {

  for (auto& block : left_->blocks()) {
    const int na = nelea_ - block.nelea;
    const int nb = neleb_ - block.neleb;

    if (na >= 0 && na <= space->norb() && nb >= 0 && nb <= space->norb()) {
      shared_ptr<RASDeterminants> det = space_->det(na, nb);
      if (det->size() > 0)
        sectors_.emplace(block, make_shared<RASBlockVectors>(det, block));
    }
  }
}


/// Copy-constructor
ProductRASCivec::ProductRASCivec(const ProductRASCivec& o) : space_(o.space_), left_(o.left_), nelea_(o.nelea_), neleb_(o.neleb_) {
  for (auto& sec : o.sectors_)
    sectors_.emplace(sec.first, make_shared<RASBlockVectors>(*sec.second));
}


/// Move-constructor
ProductRASCivec::ProductRASCivec(ProductRASCivec&& o) : sectors_(move(o.sectors_)), space_(move(o.space_)),
  left_(o.left_), nelea_(o.nelea_), neleb_(o.neleb_) {}


/// Copy-assignment
ProductRASCivec& ProductRASCivec::operator=(const ProductRASCivec& o) {
  nelea_ = o.nelea_;
  neleb_ = o.neleb_;
  left_ = o.left_;
  space_ = o.space_;
  sectors_.clear();
  for (auto& osec : o.sectors_)
    sectors_.emplace(osec.first, make_shared<RASBlockVectors>(*osec.second));
  return *this;
}


/// Move-assignment
ProductRASCivec& ProductRASCivec::operator=(ProductRASCivec&& o) {
  nelea_ = o.nelea_;
  neleb_ = o.neleb_;
  left_ = move(o.left_);
  space_ = move(o.space_);
  sectors_ = move(o.sectors_);

  return *this;
}


void ProductRASCivec::scale(const double a) {
  for (auto& sec : sectors_) sec.second->scale(a);
}


double ProductRASCivec::dot_product(const ProductRASCivec& o) const {
  assert(matches(o));
  double out = 0.0;
  for (auto& s : sectors_)
    out += s.second->dot_product(*o.sectors_.at(s.first));
  return out;
}


void ProductRASCivec::ax_plus_y(const double& a, const ProductRASCivec& o) {
  assert(matches(o));
  for (auto& s : sectors_)
    s.second->ax_plus_y(a, *o.sectors_.at(s.first));
}


void ProductRASCivec::allreduce() {
  for (auto& s : sectors_)
    s.second->allreduce();
}


void ProductRASCivec::broadcast(const int rank) {
  for (auto& s : sectors_)
    s.second->broadcast(rank);
}


void ProductRASCivec::synchronize() {
  broadcast();
}


void ProductRASCivec::print(const double thresh) const {
  for (auto& isec: sectors_) {
    bool sector_printed = false;
    const int nstate = isec.second->nstates();
    for (int ist = 0; ist < nstate; ++ist) {
      auto ci = isec.second->civec(ist);
      if (ci.norm() > thresh) {
        cout << "    " << left_->block_info_to_string(isec.first, ist) << " (x)" << endl;
        ci.print(thresh);
        sector_printed = true;
      }
    }
    if (sector_printed)
      cout << endl;
  }
}


shared_ptr<ProductRASCivec> ProductRASCivec::spin() const {
  auto out = this->clone();
  for (auto& sector : sectors()) {
    { // pure block part
      shared_ptr<Matrix> spinmat = this->left()->spin(sector.first);
      const double SAB = 0.5 * static_cast<double>((sector.first.nelea - sector.first.neleb)*sector.second->det()->nspin());
      spinmat->add_diag(SAB);
      shared_ptr<const RASBlockVectors> source = sector.second;
      shared_ptr<RASBlockVectors> target = out->sector(sector.first);
      dgemm_("N", "T", target->ndim(), target->mdim(), target->mdim(), 1.0, source->data(), source->ndim(), spinmat->data(), spinmat->ndim(),
                                                                       1.0, target->data(), target->ndim());
      // pure ras part
      const int nstate = source->mdim();
      for (int i = 0; i < nstate; ++i) {
        RAS::spin_impl(source->civec(i), target->civec(i));
      }
    }

    // mixed part
    // S_-^L S_+^RAS
    BlockKey lowered_key(sector.first.nelea-1,sector.first.neleb+1);
    if (this->contains_block(lowered_key)) {
      shared_ptr<const RASBlockVectors> source = sector.second;
      shared_ptr<RASBlockVectors> target = out->sector(lowered_key);
      shared_ptr<Matrix> spin_lower = this->left()->spin_lower(sector.first);

      BlockInfo left_state(target->left_state().nelea, target->left_state().neleb, source->mdim());
      RASBlockVectors raised_sector(target->det(), left_state);
      for (int ist = 0; ist < source->mdim(); ++ist)
        RAS::spin_raise_impl(source->civec(ist), raised_sector.civec(ist));

      dgemm_("N", "T", target->ndim(), target->mdim(), raised_sector.mdim(), 1.0, raised_sector.data(), raised_sector.ndim(), spin_lower->data(), spin_lower->ndim(),
                                                                             1.0, target->data(), target->ndim());
    }

    // S_+^L S_-^RAS
    BlockKey raised_key(sector.first.nelea+1,sector.first.neleb-1);
    if (this->contains_block(raised_key)) {
      shared_ptr<const RASBlockVectors> source = sector.second;
      shared_ptr<RASBlockVectors> target = out->sector(raised_key);
      shared_ptr<Matrix> spin_raise = this->left()->spin_raise(sector.first);

      BlockInfo left_state(target->left_state().nelea, target->left_state().neleb, source->mdim());
      RASBlockVectors lowered_sector(target->det(), left_state);
      for (int ist = 0; ist < source->mdim(); ++ist)
        RAS::spin_lower_impl(source->civec(ist), lowered_sector.civec(ist));

      dgemm_("N", "T", target->ndim(), target->mdim(), lowered_sector.mdim(), 1.0, lowered_sector.data(), lowered_sector.ndim(), spin_raise->data(), spin_raise->ndim(),
                                                                             1.0, target->data(), target->ndim());
    }
  }

  return out;
}


shared_ptr<ProductRASCivec> ProductRASCivec::spin_lower() const {
  auto out = make_shared<ProductRASCivec>(space_, left_, nelea_-1, neleb_+1);
  for (auto& source_sector : sectors()) {
    BlockKey blockkey = source_sector.first;
    BlockKey lowered_block(blockkey.nelea-1, blockkey.neleb+1);
    if (out->contains_block(lowered_block)) {
      shared_ptr<const Matrix> spin_lower_block = this->left()->spin_lower(blockkey);
      shared_ptr<Matrix> target_sector = out->sector(lowered_block);
      dgemm_("N", "T", target_sector->ndim(), target_sector->mdim(), source_sector.second->mdim(), 1.0, source_sector.second->data(), source_sector.second->ndim(),
                                                        spin_lower_block->data(), spin_lower_block->ndim(), 1.0, target_sector->data(), target_sector->ndim());
    }

    if (out->contains_block(blockkey)) {
      shared_ptr<const RASBlockVectors> target_sector = out->sector(blockkey);
      for (int ist = 0; ist < target_sector->mdim(); ++ist)
        RAS::spin_lower_impl(source_sector.second->civec(ist), target_sector->civec(ist));
    }
  }

  return out;
}


shared_ptr<ProductRASCivec> ProductRASCivec::spin_raise() const {
  auto out = make_shared<ProductRASCivec>(space_, left_, nelea_+1, neleb_-1);
  for (auto& source_sector : sectors()) {
    BlockKey blockkey = source_sector.first;
    BlockKey raised_block(blockkey.nelea+1, blockkey.neleb-1);
    if (out->contains_block(raised_block)) {
      shared_ptr<const Matrix> spin_raise_block = this->left()->spin_raise(blockkey);
      shared_ptr<Matrix> target_sector = out->sector(raised_block);
      dgemm_("N", "T", target_sector->ndim(), target_sector->mdim(), source_sector.second->mdim(), 1.0, source_sector.second->data(), source_sector.second->ndim(),
                                                        spin_raise_block->data(), spin_raise_block->ndim(), 1.0, target_sector->data(), target_sector->ndim());
    }

    if (out->contains_block(blockkey)) {
      shared_ptr<const RASBlockVectors> target_sector = out->sector(blockkey);
      for (int ist = 0; ist < target_sector->mdim(); ++ist)
        RAS::spin_raise_impl(source_sector.second->civec(ist), target_sector->civec(ist));
    }
  }

  return out;
}


double ProductRASCivec::spin_expectation() const { return this->dot_product(*spin()); }


void ProductRASCivec::spin_decontaminate(const double thresh) {
  const int nspin = nelea() - neleb();
  const int max_spin = nelea() + neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<ProductRASCivec> S2 = spin();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while( fabs(actual_expectation - pure_expectation) > thresh ) {
    if ( k > max_spin ) { this->print(0.05); throw std::runtime_error("Spin decontamination failed."); }

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    actual_expectation = dot_product(*S2);

    k += 2;
  }
}


map<BlockKey, vector<shared_ptr<ProductRASCivec>>> ProductRASCivec::split() const {
  shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(left_);
  if (!doubleblock)
    throw logic_error("ProductRASCivec::split() should only be called from a vector containing a DMRG_Block2 object.");

  shared_ptr<const DMRG_Block1> rightblock = doubleblock->right_block();
  shared_ptr<const DMRG_Block1> leftblock = doubleblock->left_block();

  set<BlockKey> used_rightblocks;
  for (auto& sec : sectors_) {
    const BlockKey& bk = sec.first;
    for (auto& bp : doubleblock->blockpairs(bk))
      used_rightblocks.insert(bp.right);
  }

  map<BlockKey, vector<shared_ptr<ProductRASCivec>>> out;
  for (BlockKey rightkey : used_rightblocks) {
    const BlockInfo rightinfo = rightblock->blockinfo(rightkey);
    vector<shared_ptr<ProductRASCivec>> outvec;

    const int nele_r  = rightinfo.nelea + rightinfo.neleb;
    const int nele_lc = (nelea_ + neleb_) - nele_r;

    for (int i = 0; i < rightinfo.nstates; ++i) {
      auto tmpout = make_shared<ProductRASCivec>(space_, leftblock, nelea_ - rightinfo.nelea, neleb_ - rightinfo.neleb);
      for (auto& outsec : tmpout->sectors()) {
        const BlockInfo leftinfo = leftblock->blockinfo(outsec.first);

        BlockKey combinedkey(rightinfo.nelea+leftinfo.nelea, rightinfo.neleb+leftinfo.neleb);
        const vector<DMRG::BlockPair>& pairs = doubleblock->blockpairs(combinedkey);
        auto iter = find_if(pairs.begin(), pairs.end(), [&leftinfo, &rightinfo] (const DMRG::BlockPair& bp)
          { return make_pair(leftinfo, rightinfo)==make_pair(bp.left, bp.right); });
        assert(iter!= pairs.end());
        const int stateindex = i*leftinfo.nstates + iter->offset;
        copy_n(sectors_.at(combinedkey)->element_ptr(0, stateindex), outsec.second->size(), outsec.second->data());

        // phase is from rearranging |CI> |L> |R> to |R> |CI> |L>
        if ((nele_lc*nele_r)%2==1)
          outsec.second->scale(-1.0);
      }
      outvec.push_back(tmpout);
    }
    out.emplace(rightkey, move(outvec));
  }

  return out;
}
