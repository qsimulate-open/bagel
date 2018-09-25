//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: subtask.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/subtask.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template class SMITH::SubTask_Merged_<2,1,5,double>;
template class SMITH::SubTask_Merged_<4,1,5,double>;
template class SMITH::SubTask_Merged_<6,1,5,double>;
template class SMITH::SubTask_Merged_<8,1,5,double>;
template class SMITH::SubTask_Merged_<2,2,5,double>;
template class SMITH::SubTask_Merged_<4,2,5,double>;
template class SMITH::SubTask_Merged_<6,2,5,double>;
template class SMITH::SubTask_Merged_<8,2,5,double>;
template class SMITH::SubTask_<1,1,double>;
template class SMITH::SubTask_<1,2,double>;
template class SMITH::SubTask_<1,3,double>;
template class SMITH::SubTask_<1,4,double>;
template class SMITH::SubTask_<1,5,double>;
template class SMITH::SubTask_<1,6,double>;
template class SMITH::SubTask_<1,7,double>;
template class SMITH::SubTask_<1,8,double>;
template class SMITH::SubTask_<1,9,double>;
template class SMITH::SubTask_<2,1,double>;
template class SMITH::SubTask_<2,2,double>;
template class SMITH::SubTask_<2,3,double>;
template class SMITH::SubTask_<2,4,double>;
template class SMITH::SubTask_<2,5,double>;
template class SMITH::SubTask_<2,6,double>;
template class SMITH::SubTask_<2,7,double>;
template class SMITH::SubTask_<2,8,double>;
template class SMITH::SubTask_<2,9,double>;
template class SMITH::SubTask_<3,1,double>;
template class SMITH::SubTask_<3,2,double>;
template class SMITH::SubTask_<3,3,double>;
template class SMITH::SubTask_<3,4,double>;
template class SMITH::SubTask_<3,5,double>;
template class SMITH::SubTask_<3,6,double>;
template class SMITH::SubTask_<3,7,double>;
template class SMITH::SubTask_<3,8,double>;
template class SMITH::SubTask_<3,9,double>;
template class SMITH::SubTask_<4,1,double>;
template class SMITH::SubTask_<4,2,double>;
template class SMITH::SubTask_<4,3,double>;
template class SMITH::SubTask_<4,4,double>;
template class SMITH::SubTask_<4,5,double>;
template class SMITH::SubTask_<4,6,double>;
template class SMITH::SubTask_<4,7,double>;
template class SMITH::SubTask_<4,8,double>;
template class SMITH::SubTask_<4,9,double>;
template class SMITH::SubTask_<5,1,double>;
template class SMITH::SubTask_<5,2,double>;
template class SMITH::SubTask_<5,3,double>;
template class SMITH::SubTask_<5,4,double>;
template class SMITH::SubTask_<5,5,double>;
template class SMITH::SubTask_<5,6,double>;
template class SMITH::SubTask_<5,7,double>;
template class SMITH::SubTask_<5,8,double>;
template class SMITH::SubTask_<5,9,double>;
template class SMITH::SubTask_<6,1,double>;
template class SMITH::SubTask_<6,2,double>;
template class SMITH::SubTask_<6,3,double>;
template class SMITH::SubTask_<6,4,double>;
template class SMITH::SubTask_<6,5,double>;
template class SMITH::SubTask_<6,6,double>;
template class SMITH::SubTask_<6,7,double>;
template class SMITH::SubTask_<6,8,double>;
template class SMITH::SubTask_<6,9,double>;
template class SMITH::SubTask_<7,1,double>;
template class SMITH::SubTask_<7,2,double>;
template class SMITH::SubTask_<7,3,double>;
template class SMITH::SubTask_<7,4,double>;
template class SMITH::SubTask_<7,5,double>;
template class SMITH::SubTask_<7,6,double>;
template class SMITH::SubTask_<7,7,double>;
template class SMITH::SubTask_<7,8,double>;
template class SMITH::SubTask_<7,9,double>;
template class SMITH::SubTask_<8,1,double>;
template class SMITH::SubTask_<8,2,double>;
template class SMITH::SubTask_<8,3,double>;
template class SMITH::SubTask_<8,4,double>;
template class SMITH::SubTask_<8,5,double>;
template class SMITH::SubTask_<8,6,double>;
template class SMITH::SubTask_<8,7,double>;
template class SMITH::SubTask_<8,8,double>;
template class SMITH::SubTask_<8,9,double>;
template class SMITH::SubTask_<9,1,double>;
template class SMITH::SubTask_<9,2,double>;
template class SMITH::SubTask_<9,3,double>;
template class SMITH::SubTask_<9,4,double>;
template class SMITH::SubTask_<9,5,double>;
template class SMITH::SubTask_<9,6,double>;
template class SMITH::SubTask_<9,7,double>;
template class SMITH::SubTask_<9,8,double>;
template class SMITH::SubTask_<9,9,double>;

template class SMITH::SubTask_<1,1,complex<double>>;
template class SMITH::SubTask_<1,2,complex<double>>;
template class SMITH::SubTask_<1,3,complex<double>>;
template class SMITH::SubTask_<1,4,complex<double>>;
template class SMITH::SubTask_<1,5,complex<double>>;
template class SMITH::SubTask_<1,6,complex<double>>;
template class SMITH::SubTask_<1,7,complex<double>>;
template class SMITH::SubTask_<1,8,complex<double>>;
template class SMITH::SubTask_<1,9,complex<double>>;
template class SMITH::SubTask_<2,1,complex<double>>;
template class SMITH::SubTask_<2,2,complex<double>>;
template class SMITH::SubTask_<2,3,complex<double>>;
template class SMITH::SubTask_<2,4,complex<double>>;
template class SMITH::SubTask_<2,5,complex<double>>;
template class SMITH::SubTask_<2,6,complex<double>>;
template class SMITH::SubTask_<2,7,complex<double>>;
template class SMITH::SubTask_<2,8,complex<double>>;
template class SMITH::SubTask_<2,9,complex<double>>;
template class SMITH::SubTask_<3,1,complex<double>>;
template class SMITH::SubTask_<3,2,complex<double>>;
template class SMITH::SubTask_<3,3,complex<double>>;
template class SMITH::SubTask_<3,4,complex<double>>;
template class SMITH::SubTask_<3,5,complex<double>>;
template class SMITH::SubTask_<3,6,complex<double>>;
template class SMITH::SubTask_<3,7,complex<double>>;
template class SMITH::SubTask_<3,8,complex<double>>;
template class SMITH::SubTask_<3,9,complex<double>>;
template class SMITH::SubTask_<4,1,complex<double>>;
template class SMITH::SubTask_<4,2,complex<double>>;
template class SMITH::SubTask_<4,3,complex<double>>;
template class SMITH::SubTask_<4,4,complex<double>>;
template class SMITH::SubTask_<4,5,complex<double>>;
template class SMITH::SubTask_<4,6,complex<double>>;
template class SMITH::SubTask_<4,7,complex<double>>;
template class SMITH::SubTask_<4,8,complex<double>>;
template class SMITH::SubTask_<4,9,complex<double>>;
template class SMITH::SubTask_<5,1,complex<double>>;
template class SMITH::SubTask_<5,2,complex<double>>;
template class SMITH::SubTask_<5,3,complex<double>>;
template class SMITH::SubTask_<5,4,complex<double>>;
template class SMITH::SubTask_<5,5,complex<double>>;
template class SMITH::SubTask_<5,6,complex<double>>;
template class SMITH::SubTask_<5,7,complex<double>>;
template class SMITH::SubTask_<5,8,complex<double>>;
template class SMITH::SubTask_<5,9,complex<double>>;
template class SMITH::SubTask_<6,1,complex<double>>;
template class SMITH::SubTask_<6,2,complex<double>>;
template class SMITH::SubTask_<6,3,complex<double>>;
template class SMITH::SubTask_<6,4,complex<double>>;
template class SMITH::SubTask_<6,5,complex<double>>;
template class SMITH::SubTask_<6,6,complex<double>>;
template class SMITH::SubTask_<6,7,complex<double>>;
template class SMITH::SubTask_<6,8,complex<double>>;
template class SMITH::SubTask_<6,9,complex<double>>;
template class SMITH::SubTask_<7,1,complex<double>>;
template class SMITH::SubTask_<7,2,complex<double>>;
template class SMITH::SubTask_<7,3,complex<double>>;
template class SMITH::SubTask_<7,4,complex<double>>;
template class SMITH::SubTask_<7,5,complex<double>>;
template class SMITH::SubTask_<7,6,complex<double>>;
template class SMITH::SubTask_<7,7,complex<double>>;
template class SMITH::SubTask_<7,8,complex<double>>;
template class SMITH::SubTask_<7,9,complex<double>>;
template class SMITH::SubTask_<8,1,complex<double>>;
template class SMITH::SubTask_<8,2,complex<double>>;
template class SMITH::SubTask_<8,3,complex<double>>;
template class SMITH::SubTask_<8,4,complex<double>>;
template class SMITH::SubTask_<8,5,complex<double>>;
template class SMITH::SubTask_<8,6,complex<double>>;
template class SMITH::SubTask_<8,7,complex<double>>;
template class SMITH::SubTask_<8,8,complex<double>>;
template class SMITH::SubTask_<8,9,complex<double>>;
template class SMITH::SubTask_<9,1,complex<double>>;
template class SMITH::SubTask_<9,2,complex<double>>;
template class SMITH::SubTask_<9,3,complex<double>>;
template class SMITH::SubTask_<9,4,complex<double>>;
template class SMITH::SubTask_<9,5,complex<double>>;
template class SMITH::SubTask_<9,6,complex<double>>;
template class SMITH::SubTask_<9,7,complex<double>>;
template class SMITH::SubTask_<9,8,complex<double>>;
template class SMITH::SubTask_<9,9,complex<double>>;

#endif
