//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relgrad_base.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#ifndef __SRC_MAT1E_GRAD_RELGRAD_BASE_H
#define __SRC_MAT1E_GRAD_RELGRAD_BASE_H

#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>
#include <src/util/math/matrix.h>
#include <src/util/math/xyzfile.h>
#include <src/wfn/geometry.h>

namespace bagel {

	class Relgrad_base {

	protected:
		int natom;
		int nbasis;
	  	int nunc;

	    std::shared_ptr<const Geometry> mol;
	    std::shared_ptr<const Molecule> molu;

	    const MixedBasis<OverlapBatch> U_T;
	    std::vector<std::vector<Matrix>> PU;

	    std::vector<std::vector<Matrix>> S_X;
	    std::vector<std::vector<Matrix>> T_X;
	    std::vector<std::vector<Matrix>> V_X;
	    std::vector<std::vector<Matrix>> O_X;

	    const Matrix id;
	    std::map<std::shared_ptr<VectorB>, std::shared_ptr<Matrix>> vec2mat;

	    void store_mat(const shared_ptr<VectorB>);
	    void contracts(const shared_ptr<Matrix>);
	    void overlapgrad(const shared_ptr<Matrix>);
	    void kineticgrad(const shared_ptr<Matrix>);
	    void naigrad(const shared_ptr<Matrix>);
	    void smallnaigrad(const shared_ptr<Matrix>);

	public:
	    Relgrad_base() { }
	    Relgrad_base(std::shared_ptr<const Geometry>);

	};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Relgrad_base)

#endif

