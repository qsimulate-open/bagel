//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/rel/dfhalfcomplex.h>

using namespace std;
using namespace bagel;

DFHalfComplex::DFHalfComplex(const shared_ptr<const DFDist> df, shared_ptr<const Matrix> rcoeff, shared_ptr<const Matrix> icoeff, 
                             const bool swap, pair<const int, const int> coord, pair<const int, const int> basis) : coord_(coord), basis_(basis) {

  dim_ = rcoeff->ndim();
  shared_ptr<DFHalfDist> rhalfbj;
  shared_ptr<DFHalfDist> ihalfbj;

  if (swap == true) {
    rhalfbj = df->compute_half_transform_swap(rcoeff);
    ihalfbj = df->compute_half_transform_swap(icoeff); 
  } else {
    rhalfbj = df->compute_half_transform(rcoeff);
    ihalfbj = df->compute_half_transform(icoeff); 
  }

  dfhalf_[0] = rhalfbj->apply_J();
  dfhalf_[1] = ihalfbj->apply_J();
}

#if 0
void DFHalfComplex::add_Exop_block(shared_ptr<ZMatrix> dfock, shared_ptr<DFHalfComplex> dfc) {

  complex<double> coeff1 = compute_coeff(dfc->basis(), dfc->coord());
  double coeff2 = (coeff1.real() == 0.0 ? -1.0 :1.0);
  const tuple<int, int, int, int> index = compute_index_Exop(dfc->basis(), dfc->coord());

  array<shared_ptr<Matrix>, 4> exop;
  *exop[0] = *dfhalf_[0]->form_2index(dfc->get_real(), -1.0); 
  *exop[1] = *dfhalf_[1]->form_2index(dfc->get_imag(), -1.0);
  *exop[2] = *dfhalf_[0]->form_2index(dfc->get_imag(), coeff2); 
  *exop[3] = *dfhalf_[1]->form_2index(dfc->get_real(), -coeff2);

  for (int i = 0; i != 4; ++i) {
    double factor = (i == 2 ? -1.0 : 1.0);
    shared_ptr<Matrix> texop = exop[i]->transpose(factor);
    dfock->add_real_block(coeff1, dim_ * get<0>(index), dim_ * get<1>(index), dim_, dim_, exop[i]);
    if (get<2>(index) != -1)
      dfock->add_real_block(coeff1, dim_ * get<2>(index), dim_ * get<3>(index), dim_, dim_, texop);
    //cross terms
    if (dfc->basis().second != basis_.second) {
      dfock->add_real_block(coeff1, dim_ * get<1>(index), dim_ * get<0>(index), dim_, dim_, texop);
      if (get<2>(index) != -1)
        dfock->add_real_block(coeff1, dim_ * get<3>(index), dim_ * get<2>(index), dim_, dim_, exop[i]);
    }
  }
}

void DFHalfComplex::add_Jop_block(shared_ptr<ZMatrix> dfock, shared_ptr<const DFData> dfdata, shared_ptr<const Matrix> trocoeff, 
                 shared_ptr<const Matrix> tiocoeff) {

  shared_ptr<const DFDist> df = dfdata->df();
  
  complex<double> coeff1 = compute_coeff(dfdata->basis(), dfdata->coord());
  double coeff2 = (coeff1.real() == 0.0 ? -1.0 :1.0);
  const tuple<int, int, int, int> index = compute_index_Jop(dfdata->basis(), dfdata->coord());

  array<shared_ptr<Matrix>, 4> jop;
  *jop[0] = *df->compute_Jop(dfhalf_[0], trocoeff, true);
  *jop[1] = *df->compute_Jop(dfhalf_[1], tiocoeff, true);
  *jop[2] = *df->compute_Jop(dfhalf_[0], tiocoeff, true) * (-1.0) * coeff2;
  *jop[3] = *df->compute_Jop(dfhalf_[1], trocoeff, true) * coeff2;

  for (int i = 0; i != 4; ++i) {
   //add it twice, once to first basis combo, then once to opposite basis combo
    dfock->add_real_block(coeff1, dim_ * get<0>(index), dim_ * get<1>(index), dim_, dim_, jop[i]);
    dfock->add_real_block(coeff1, dim_ * get<2>(index), dim_ * get<3>(index), dim_, dim_, jop[i]);

    //if basis1 != basis2, get transpose to fill in opposite corner
    if (dfdata->cross()) {
      shared_ptr<Matrix> tjop = jop[i]->transpose();
      dfock->add_real_block(coeff1, dim_ * get<1>(index), dim_ * get<0>(index), dim_, dim_, (*tjop * dfdata->cross_coeff()).data());
      dfock->add_real_block(coeff1, dim_ * get<3>(index), dim_ * get<2>(index), dim_, dim_, (*tjop * dfdata->cross_coeff()).data());
    }
  }
}
#endif

complex<double> DFHalfComplex::compute_coeff(pair<const int, const int> basis2, pair<const int, const int> coord2) {
  const double tc = 1.0 / (2.0* c__);
  double power = 0.0;
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);

  // Xa Xb, Ya Yb, Za Zb respectively
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};

  complex<double> coeff1 = rcoeff, coeff2 = rcoeff, coeff3 = rcoeff, coeff4 = rcoeff;

  // |coord_ coord_*)(coord2* coord2)
  if (coord_.first != -1) {
    power+=2.0;
    coeff1 = coeff[coord_.first][basis_.first];
    coeff2 = (coord_.second == 1) ? -coeff[coord_.second][basis_.second] : coeff[coord_.second][basis_.second]; 
  }
  if (coord2.first != -1) {
    power+=2.0;
    coeff3 = (coord2.first == 1) ? -coeff[coord2.first][basis2.first] : coeff[coord2.first][basis2.first];
    coeff4 = coeff[coord2.second][basis2.second];
  }

  complex<double> out = coeff1 * coeff2 * coeff3 * coeff4 * ::pow(tc, power);
  return out;
}

const tuple<int, int, int, int> DFHalfComplex::compute_index_Exop(pair<const int, const int> basis2, pair<const int, const int> coord2) {

  int index1, index2, index3, index4;

  // 4x4 ZMatrix starting at 0,0 (large, large) or 0,2n (large, small) or 2n,0 (small, large) or 2n,2n (small)
  int start1 = (coord_.first == -1 ? 0 : 2);
  int start2 = (coord2.first == -1 ? 0 : 2);
  //go from small large to large small or vice versa
  index1 = start1 + basis_.second;
  index2 = start2 + basis2.second;

  if (start1 != start2) {
#if 0
    int opp1 = (coord_.first == -1 ? 2 : 0);
    int opp2 = (coord2.first == -1 ? 2 : 0);
    index3 = opp1 + basis_.second;
    index4 = opp2 + basis2.second;
#else
    index3 = index2;
    index4 = index1;
#endif
  } else {
    index3 = -1;
    index4 = -1;
  }
  
  return make_tuple(index1, index2, index3, index4);
}

const tuple<int, int, int, int> DFHalfComplex::compute_index_Jop(pair<const int, const int> basis, pair<const int, const int> coord) {
  int index1, index2, index3, index4;
  // 4x4 ZMatrix either starting at 0,0 (large) or 2n,2n (small)
  int start = (coord.first == -1 ? 0 : 2);
  // put transposed Matrices in submatrix opposite original
  int opp1 =  (basis.first  == 0 ? 1 : 0);
  int opp2 =  (basis.second == 0 ? 1 : 0);

  index1 = start + basis.first;
  index2 = start + basis.second;
  index3 = start + opp1;
  index4 = start + opp2; 

  return make_tuple(index1, index2, index3, index4);
}

const int DFHalfComplex::coeff_matrix() const {
  
  int out = (coord_.first == -1 ? basis_.second : basis_.second + 2);

  return out;
}
