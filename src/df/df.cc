//
// Newint - Parallel electron correlation program.
// Filename: df.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <memory>
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/f77.h>
#include <src/df/df.h>
#include <src/rysint/eribatch.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <list>

using namespace std;


class DFIntTask {
  protected:
    vector<shared_ptr<const Shell> > shell_;
    vector<int> offset_;
    DensityFit* df_;

  public:
    DFIntTask(vector<shared_ptr<const Shell> > a, vector<int> b, DensityFit* df) : shell_(a), offset_(b), df_(df) {};    
    DFIntTask() {};

    void compute() {
      const double* ppt = df_->compute_batch(shell_);

      const size_t nbasis1 = df_->nbasis1();
      const size_t naux = df_->naux();
      // all slot in
      if (offset_.size() == 3) {
        double* data = df_->data();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {  
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1) {  
            for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++ppt) {  
              data[j2+naux*(j1+nbasis1*j0)] = data[j2+naux*(j0+nbasis1*j1)] = *ppt;
            }
          }
        }
      } else if (offset_.size() == 2) {
        double* data = df_->data2();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0) {  
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {  
            data[j1+j0*naux] = data[j0+j1*naux] = *ppt;
          }
        }
      } else {
        assert(false);
      }
    };
};


void DensityFit::common_init(const vector<shared_ptr<const Atom> >& atoms0,  const vector<vector<int> >& offsets0,
                             const vector<shared_ptr<const Atom> >& atoms1,  const vector<vector<int> >& offsets1,
                             const vector<shared_ptr<const Atom> >& aux_atoms,  const vector<vector<int> >& aux_offsets, const double throverlap, 
                             const bool compute_inverse) {

  // this will be distributed in the future.
  data_ = unique_ptr<double[]>(new double[nbasis0_*nbasis1_*naux_]);
  fill(data_.get(), data_.get()+nbasis0_*nbasis1_*naux_, 0.0);

  const shared_ptr<const Shell> b3(new Shell(atoms0.front()->shells().front()->spherical()));

  if (atoms0 == atoms1) {
    // these atom loops will be distributed across nodes
    // TODO align so that heavy atoms come first

    auto oa0 = offsets0.begin();
    for (auto a0 = atoms0.begin(); a0 != atoms0.end(); ++a0, ++oa0) {
      auto oa1 = oa0;
      for (auto a1 = a0; a1 != atoms0.end(); ++a1, ++oa1) {
        auto oa2 = aux_offsets.begin(); 
        for (auto a2 = aux_atoms.begin(); a2 != aux_atoms.end(); ++a2, ++oa2) {

          // make a list of input and offsets
          list<DFIntTask> tasks;
          auto o0 = oa0->begin();
          for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
            auto o1 = a0!=a1 ? oa1->begin() : o0;
            for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {
              auto o2 = oa2->begin();
              for (auto b2 = (*a2)->shells().begin(); b2 != (*a2)->shells().end(); ++b2, ++o2) {
                tasks.push_back(DFIntTask(vector<shared_ptr<const Shell> >{{b3, *b2, *b1, *b0}}, vector<int>{{*o0, *o1, *o2}}, this));
              }
            }
          }
                
          // these shell loops will be distributed across threads 
          const int num_threads = 2;
          TaskQueue<DFIntTask> tq(tasks);
          tq.compute(num_threads);
        }
      }
    }
  } else {
    throw logic_error("Dual basis DF builder has not been implemented yet. See src/df/df.cc");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  data2_ = unique_ptr<double[]>(new double[naux_*naux_]);
  fill(data2(), data2()+naux_*naux_, 0.0);

  // will be distributed across nodes
  auto oa0 = aux_offsets.begin();
  for (auto a0 = aux_atoms.begin(); a0 != aux_atoms.end(); ++a0, ++oa0) {
    auto oa1 = oa0;
    for (auto a1 = a0; a1 != aux_atoms.end(); ++a1, ++oa1) {

      // make a list of input and offsets
      list<DFIntTask> tasks;
      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = a0!=a1 ? oa1->begin() : o0;
        for (auto b1 = (a0!=a1 ? (*a1)->shells().begin() : b0); b1 != (*a1)->shells().end(); ++b1, ++o1) {
          tasks.push_back(DFIntTask(vector<shared_ptr<const Shell> >{{*b1, b3, *b0, b3}}, vector<int>{{*o0, *o1}}, this));
        }
      }

      // these shell loops will be distributed across threads 
      for (auto itask = tasks.begin(); itask != tasks.end(); ++itask)
        itask->compute(); 

    }
  }

  if (compute_inverse) {
    const size_t lwork = 5*naux_;
    unique_ptr<double[]> vec(new double[naux_]);
    unique_ptr<double[]> work(new double[max(lwork,naux_*naux_)]);

    int info;
    dsyev_("V", "U", naux_, data2_, naux_, vec, work, lwork, info); 
    if (info) throw runtime_error("dsyev failed in DF fock builder");

    for (int i = 0; i != naux_; ++i)
      vec[i] = vec[i] > throverlap ? 1.0/sqrt(sqrt(vec[i])) : 0.0;
    for (int i = 0; i != naux_; ++i) {
      for (int j = 0; j != naux_; ++j) {
        data2_[j+i*naux_] *= vec[i];
      }
    }
    // work now contains -1/2
    dgemm_("N", "T", naux_, naux_, naux_, 1.0, data2_, naux_, data2_, naux_, 0.0, work, naux_); 
    copy(work.get(), work.get()+naux_*naux_, data2());
  }

}


unique_ptr<double[]> DensityFit::compute_Jop(const double* den) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  unique_ptr<double[]> tmp0 = compute_cd(den);
  unique_ptr<double[]> out(new double[nbasis0_*nbasis1_]);
  // then compute J operator J_{rs} = |E*) (E|rs)
  dgemv_("T", naux_, nbasis0_*nbasis1_, 1.0, data_, naux_, tmp0, 1, 0.0, out, 1); 
  return out;
} 


unique_ptr<double[]> DensityFit::compute_cd(const double* den) const {
  unique_ptr<double[]> tmp0(new double[naux_]);
  unique_ptr<double[]> tmp1(new double[naux_]);
  dgemv_("N", naux_, nbasis0_*nbasis1_, 1.0, data_.get(), naux_, den, 1, 0.0, tmp0.get(), 1); 
  dgemv_("N", naux_, naux_, 1.0, data2_, naux_, tmp0, 1, 0.0, tmp1, 1); 
  dgemv_("N", naux_, naux_, 1.0, data2_, naux_, tmp1, 1, 0.0, tmp0, 1); 
  return tmp0;
}


DF_AO::DF_AO(const int nbas0, const int nbas1, const int naux, const vector<const double*> cd, const vector<const double*> dd)
 : DensityFit(nbas0, nbas1, naux) {
  assert(cd.size() == dd.size());

  // initialize to zero
  unique_ptr<double[]> buf(new double[naux_*nbasis0_*nbasis1_]); 
  fill(buf.get(), buf.get()+size(), 0.0);

  for (auto citer = cd.begin(), diter = dd.begin(); citer != cd.end(); ++citer, ++diter) {
    dger_(naux_, nbasis0_*nbasis1_, 1.0, *citer, 1, *diter, 1, buf.get(), naux_);
  }
  data_ = move(buf);
}


shared_ptr<DF_Half> DF_Half::clone() const {
  unique_ptr<double[]> dat(new double[size()]);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, dat));
}


shared_ptr<DF_Half> DF_Half::copy() const {
  unique_ptr<double[]> dat(new double[size()]);
  std::copy(data_.get(), data_.get()+size(), dat.get());
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, dat));
}


shared_ptr<DF_Half> DF_Half::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
  dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, d->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, out));
}


shared_ptr<DF_Half> DF_Half::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Half)");
  unique_ptr<double[]> jj(new double[naux_*naux_]);
  dgemm_("N", "N", naux_, naux_, naux_, 1.0, d->data_2index(), naux_, d->data_2index(), naux_, 0.0, jj.get(), naux_);

  unique_ptr<double[]> out(new double[nocc_*naux_*nbasis_]);
  dgemm_("N", "N", naux_, nocc_*nbasis_, naux_, 1.0, jj, naux_, data_, naux_, 0.0, out, naux_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, out));
}


shared_ptr<DF_Half> DensityFit::compute_half_transform(const double* c, const size_t nocc) const {
  // it starts with DGEMM of inner index (for some reasons)...
  unique_ptr<double[]> tmp(new double[naux_*nbasis0_*nocc]);
  for (size_t i = 0; i != nbasis0_; ++i) {
    dgemm_("N", "N", naux_, nocc, nbasis1_, 1.0, data_.get()+i*naux_*nbasis1_, naux_, c, nbasis1_, 0.0, tmp.get()+i*naux_*nocc, naux_);
  }
  return shared_ptr<DF_Half>(new DF_Half(shared_from_this(), nocc, tmp));
}


void DF_Half::form_2index(unique_ptr<double[]>& target, const double a, const double b) const {
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, nbasis_, common, a, data_.get(), common, data_.get(), common, b, target.get(), nbasis_); 
}


unique_ptr<double[]> DF_Half::form_2index(shared_ptr<const DF_Full> o, const double a, const double b) const {
  assert(b == 0);
  unique_ptr<double[]> tmp(new double[nbasis_*o->nocc2()]);
  form_2index(tmp, o, a, b);
  return tmp;
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DF_Full> o, const double a, const double b) const {
  if (nocc_ != o->nocc1()) throw logic_error("nocc_ and o->nocc1() should be the same: in DF_Half::form_2index");
  const int common = nocc_ * naux_;
  dgemm_("T", "N", nbasis_, o->nocc2(), common, a, data_.get(), common, o->data(), common, b, target.get(), nbasis_); 
}


void DF_Half::form_2index(unique_ptr<double[]>& target, shared_ptr<const DensityFit> o, const double a) const {
  fill(target.get(), target.get()+nocc_*nbasis_, 0.0);
  for (size_t i = 0; i != nbasis_; ++i)
    dgemm_("T", "N", nocc_, nbasis_, naux_, a, data_.get()+i*naux_*nocc_, naux_, o->data_3index()+i*naux_*nbasis_, naux_,
                                            1.0, target.get(), nocc_);
}
unique_ptr<double[]> DF_Half::form_2index(shared_ptr<const DensityFit> o, const double a) const {
  unique_ptr<double[]> out(new double[nocc_*nbasis_]);
  form_2index(out, o, a);
  return move(out);
}


void DF_Half::form_4index(unique_ptr<double[]>& target) const {
  const int ndim = nbasis_ * nocc_;
  dgemm_("T", "N", ndim, ndim, naux_, 1.0, data_.get(), naux_, data_.get(), naux_, 0.0, target.get(), ndim); 
}

unique_ptr<double[]> DF_Half::form_4index() const {
  const size_t ndim = nbasis_ * nocc_;
  unique_ptr<double[]> out(new double[ndim*ndim]);
  form_4index(out);
  return move(out);
}

shared_ptr<DF_Full> DF_Half::compute_second_transform(const double* c, const size_t nocc) const {
  unique_ptr<double[]> tmp(new double[naux_*nocc_*nocc]);
  dgemm_("N", "N", naux_*nocc_, nocc, nbasis_, 1.0, data_.get(), naux_*nocc_, c, nbasis_, 0.0, tmp.get(), naux_*nocc_); 
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc_, nocc, tmp)); 
}


unique_ptr<double[]> DF_Half::compute_Kop_1occ(const double* den) const {
  return apply_density(den)->form_2index(df_);
}


unique_ptr<double[]> DF_Half::form_aux_2index(const shared_ptr<const DF_Half> o) const {
  unique_ptr<double[]> out(new double[naux_*naux_]);
  if (nocc_*nbasis_ != o->nocc_*o->nbasis_) throw logic_error("wrong call to DF_Full::form_aux_2index");
  dgemm_("N", "T", naux_, naux_, nocc_*nbasis_, 1.0, data(), naux_, o->data(), naux_, 0.0, out.get(), naux_); 
  return out;
}



shared_ptr<DF_Half> DF_Half::apply_density(const double* den) const {
  unique_ptr<double[]> buf(new double[naux_*nbasis_*nocc_]);
  dgemm_("N", "N", naux_*nocc_, nbasis_, nbasis_, 1.0, data(), naux_*nocc_, den, nbasis_, 0.0, buf.get(), naux_*nocc_);
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc_, buf));
}


shared_ptr<DF_Full> DF_Full::apply_J(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, d->data_2index(), naux_, data_.get(), naux_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_JJ(shared_ptr<const DensityFit> d) const {
  if (!d->data_2index()) throw logic_error("apply_J called from an object without a 2 index integral (DF_Full)");
  unique_ptr<double[]> jj(new double[naux_*naux_]);
  dgemm_("N", "N", naux_, naux_, naux_, 1.0, d->data_2index(), naux_, d->data_2index(), naux_, 0.0, jj.get(), naux_);

  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "N", naux_, nocc1_*nocc2_, naux_, 1.0, jj, naux_, data_, naux_, 0.0, out, naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm) const {
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  dgemm_("N", "T", naux_, nocc1_*nocc2_, nocc1_*nocc2_, 1.0, data_.get(), naux_, rdm, nocc1_*nocc2_, 0.0, out.get(), naux_);
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


shared_ptr<DF_Full> DF_Full::apply_2rdm(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == nocc1_ && nocc1_ == nocc2_);
  // checking if natural orbitals...
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__) throw logic_error("DF_Full::apply_2rdm should be called with natural orbitals");
  }
  unique_ptr<double[]> out(new double[nocc1_*nocc2_*naux_]);
  fill(out.get(), out.get()+nocc1_*nocc2_*naux_, 0.0);
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      daxpy_(naux_, -2.0, data_.get()+naux_*(j+nocc1_*i), 1, out.get()+naux_*(j+nocc1_*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 1.0, data()+naux_*(i+nocc1_*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 4.0, diagsum.get(), 1, out.get()+naux_*(i+nocc1_*i), 1);

  // act-act part
  // compress
  unique_ptr<double[]> buf(new double[nact*nact*naux_]);
  unique_ptr<double[]> buf2(new double[nact*nact*naux_]);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(naux_, data_.get()+naux_*(j+nclosed+nocc1_*(i+nclosed)), 1, buf.get()+naux_*(j+nact*i),1); 
  // multiply
  dgemm_("N", "N", naux_, nact*nact, nact*nact, 1.0, buf.get(), naux_, rdm, nact*nact, 0.0, buf2.get(), naux_);
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(naux_, buf2.get()+naux_*(j+nact*i),1, out.get()+naux_*(j+nclosed+nocc1_*(i+nclosed)), 1);

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab 
  // ASSUMING natural orbitals
  for (int i = 0; i != nact; ++i)
    daxpy_(naux_, 2.0*rdm1[i+nact*i], diagsum.get(), 1, out.get()+naux_*(i+nclosed+nocc1_*(i+nclosed)), 1); 
  unique_ptr<double[]> diagsum2(new double[naux_]);
  dgemv_("N", naux_, nact*nact, 1.0, buf.get(), naux_, rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(naux_, 2.0, diagsum2.get(), 1, out.get()+naux_*(i+nocc1_*i), 1); 
  // exchange contribution
  for (int i = 0; i != nact; ++i) { 
    for (int j = 0; j != nclosed; ++j) { 
      daxpy_(naux_, -rdm1[i+nact*i], data_.get()+naux_*(j+nocc1_*(i+nclosed)), 1, out.get()+naux_*(j+nocc1_*(i+nclosed)), 1);
      daxpy_(naux_, -rdm1[i+nact*i], data_.get()+naux_*(i+nclosed+nocc1_*j), 1, out.get()+naux_*(i+nclosed+nocc1_*j), 1);
    }
  }

  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, out));
}


// forms all-internal 4-index MO integrals
void DF_Full::form_4index(unique_ptr<double[]>& target) const {
  const int dim = nocc1_ * nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, dim, naux, 1.0, data_.get(), naux, data_.get(), naux, 0.0, target.get(), dim); 
}


unique_ptr<double[]> DF_Full::form_4index() const {
  unique_ptr<double[]> out(new double[nocc1_ * nocc2_ * nocc1_ * nocc2_]);
  form_4index(out);
  return move(out);
}


void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nocc1_ * o->nocc2_;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, odim, naux, 1.0, data_.get(), naux, o->data_.get(), naux, 0.0, target.get(), dim); 
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DF_Full> o) const {
  unique_ptr<double[]> out(new double[nocc1_ * nocc2_ * o->nocc1_ * o->nocc2_]);
  form_4index(out, o);
  return move(out);
}


// for MP2-like quantities
void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o, const size_t n) const {
  const int dim = nocc1_ * nocc2_;
  const int odim = o->nocc1_; // o->nocc2_ is fixed at n;
  const int naux = df_->naux();
  dgemm_("T", "N", dim, odim, naux, 1.0, data_.get(), naux, o->data_.get()+naux*odim*n, naux, 0.0, target.get(), dim); 
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DF_Full> o, const size_t n) const {
  unique_ptr<double[]> target(new double[o->nocc1_*nocc1_*nocc2_]);
  form_4index(target, o, n);
  return move(target);
}


// Joperator. Note that (r,s) runs first; i.e., in the operator form
void DF_Full::form_4index(unique_ptr<double[]>& target, const shared_ptr<const DensityFit> o) const {
  const size_t dim = nocc1_ * nocc2_;
  const size_t odim = o->nbasis0() * o->nbasis1();
  shared_ptr<DF_Full> tmp = this->apply_J();
  dgemm_("T", "N", odim, dim, naux_, 1.0, o->data_3index(), naux_, tmp->data(), naux_, 0.0, target.get(), odim); 
}

unique_ptr<double[]> DF_Full::form_4index(const shared_ptr<const DensityFit> o) const {
  const size_t dim0 = nocc1_ * o->nbasis0();
  const size_t dim1 = nocc2_ * o->nbasis1();
  assert(nocc1_ == nocc2_); // TODO <- do you need this?
  unique_ptr<double[]> out(new double[dim0*dim1]);
  form_4index(out, o);
  return move(out);
}


unique_ptr<double[]> DF_Full::form_aux_2index(const shared_ptr<const DF_Full> o) const {
  unique_ptr<double[]> out(new double[naux_*naux_]);
  if (nocc1_*nocc2_ != o->nocc1_*o->nocc2_) throw logic_error("wrong call to DF_Full::form_aux_2index");
  dgemm_("N", "T", naux_, naux_, nocc1_*nocc2_, 1.0, data(), naux_, o->data(), naux_, 0.0, out.get(), naux_); 
  return out;
}


void DF_Full::form_2index(unique_ptr<double[]>& target, const shared_ptr<const DF_Half> o, const double a) {
  assert(nocc1() == o->nocc());
  dgemm_("T", "N", nocc2_, o->nbasis(), nocc1_*naux_, a, data(), nocc1_*naux_, o->data(), nocc1_*naux_, 0.0, target.get(), nocc2_); 
}


unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Half> o, const double a) {
  unique_ptr<double[]> out(new double[nocc2_*o->nbasis()]);
  form_2index(out, o, a);
  return out;
}


void DF_Full::form_2index(unique_ptr<double[]>& target, const shared_ptr<const DF_Full> o, const double a) {
  assert(nocc1_ == o->nocc1_);
  dgemm_("T", "N", nocc2_, o->nocc2(), nocc1_*naux_, a, data(), nocc1_*naux_, o->data(), nocc1_*naux_, 0.0, target.get(), nocc2_); 
}


unique_ptr<double[]> DF_Full::form_2index(const shared_ptr<const DF_Full> o, const double a) {
  unique_ptr<double[]> out(new double[nocc2_*o->nocc2()]);
  form_2index(out, o, a);
  return out;
}


shared_ptr<DF_Full> DF_Full::clone() const {
  unique_ptr<double[]> d(new double[size()]);
  std::fill(d.get(), d.get()+size(), 0.0); 
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


shared_ptr<DF_Full> DF_Full::copy() const {
  unique_ptr<double[]> d(new double[size()]);
  std::copy(data(), data()+size(), d.get()); 
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


void DF_Full::daxpy(const double a, const DF_Full& o) {
  daxpy_(size(), a, o.data(), 1, data(), 1);
}


void DF_Full::daxpy(const double a, const DF_Half& o) {
  if (o.size() != size()) throw logic_error("DF_Full::daxpy was called in a wrong way...");
  daxpy_(size(), a, o.data(), 1, data(), 1);
}


void DF_Full::scale(const double a) {
  dscal_(size(), a, data(), 1);
}


// TODO THIS FUNCTION IS VERY INEFFICIENT
// note that this function symmetrizes but not divides by 2
void DF_Full::symmetrize() {
  assert(nocc1_ == nocc2_);
  const int n = nocc1_;
  for (int i = 0; i != n; ++i) {
    for (int j = i; j != n; ++j) {
      for (int k = 0; k != naux_; ++k) {
        data_[k+naux_*(j+n*i)] = data_[k+naux_*(i+n*j)] = (data_[k+naux_*(j+n*i)] + data_[k+naux_*(i+n*j)]);
      }
    }
  } 
}


shared_ptr<DF_Full> DF_Full::apply_closed_2RDM() const {
  assert(nocc1_ == nocc2_);
  const int nocc = nocc1_; 
  unique_ptr<double[]> d(new double[size()]);
  fill(d.get(), d.get()+size(), 0.0);
  // exchange contributions
  daxpy_(size(), -2.0, data(), 1, d.get(), 1);
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, 1.0, data()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
  }
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, 4.0, diagsum.get(), 1, d.get()+naux_*(i+nocc*i), 1);
  }
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
}


// Caution
//   o strictly assuming that we are using natural orbitals. 
//
shared_ptr<DF_Full> DF_Full::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(nocc1_ == nocc2_);
  const int nocc = nocc1_; 
  unique_ptr<double[]> d(new double[size()]);
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data(), naux_*nocc, amat, nocc, 0.0, d2.get(), naux_*nocc);
    for (int i = 0; i != nocc; ++i) {
      dgemm_("N", "N", naux_, nocc, nocc, -1.0, d2.get()+naux_*nocc*i, naux_, amat, nocc, 0.0, d.get()+naux_*nocc*i, naux_);
    }
    dgemm_("N", "N", naux_*nocc, nocc, nocc, 1.0, data(), naux_*nocc, bmat, nocc, 0.0, d2.get(), naux_*nocc);
    for (int i = 0; i != nocc; ++i) {
      dgemm_("N", "N", naux_, nocc, nocc, -1.0, d2.get()+naux_*nocc*i, naux_, bmat, nocc, 1.0, d.get()+naux_*nocc*i, naux_);
    }
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc]; 
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[naux_]);
  fill(diagsum.get(), diagsum.get()+naux_, 0.0);
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, sum[i], data()+naux_*(i+nocc*i), 1, diagsum.get(), 1);
  }
  for (int i = 0; i != nocc; ++i) {
    daxpy_(naux_, sum[i], diagsum.get(), 1, d.get()+naux_*(i+nocc*i), 1);
  }
  return shared_ptr<DF_Full>(new DF_Full(df_, nocc1_, nocc2_, d));
} 


// AO back transformation (q|rs)[CCdag]_rt [CCdag]_su 
shared_ptr<DF_Half> DF_Full::back_transform(const double* c) const{
  const int nbas = df_->nbasis1();
  unique_ptr<double[]> d(new double[nbas*nocc1_*naux_]);
  dgemm_("N", "T", naux_*nocc1_, nbas, nocc2_, 1.0, data(), naux_*nocc1_, c, nbas, 0.0, d.get(), naux_*nocc1_); 
  return shared_ptr<DF_Half>(new DF_Half(df_, nocc1_, d)); 
}


shared_ptr<DF_AO> DF_Half::back_transform(const double* c) const{
  const int nbas = df_->nbasis0();
  unique_ptr<double[]> d(new double[nbas*nbasis_*naux_]);
  for (int i = 0; i != nbasis_; ++i)
    dgemm_("N", "T", naux_, nbas, nocc_, 1.0, data()+i*naux_*nocc_, naux_, c, nbas, 0.0, d.get()+i*naux_*nbas, naux_); 
  return shared_ptr<DF_AO>(new DF_AO(nbas, nbasis_, naux_, d)); 
}


