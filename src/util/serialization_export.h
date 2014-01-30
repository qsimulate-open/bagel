//
// BAGEL - Parallel electron correlation program.
// Filename: serialization_export.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

//
// **** CAUTION ****
// this file should be included once!
//

#ifdef SERIALIZATION_EXPORT_IN_MAIN

#include <boost/serialization/export.hpp>
#include <src/scf/scf.h>
#include <src/scf/uhf.h>
#include <src/scf/rohf.h>
#include <src/scf/fock.h>
#include <src/molecule/kinetic.h>
#include <src/molecule/hcore.h>
#include <src/molecule/overlap.h>
#include <src/math/matrix.h>
#include <src/math/zmatrix.h>
#include <src/math/diis.h>

// SCF
BOOST_CLASS_EXPORT(bagel::SCF);
BOOST_CLASS_EXPORT(bagel::UHF);
BOOST_CLASS_EXPORT(bagel::ROHF);
BOOST_CLASS_EXPORT(bagel::Fock<0>);
BOOST_CLASS_EXPORT(bagel::Fock<1>);

// Molecule
BOOST_CLASS_EXPORT(bagel::Overlap);
BOOST_CLASS_EXPORT(bagel::Hcore);
BOOST_CLASS_EXPORT(bagel::Matrix1e);
BOOST_CLASS_EXPORT(bagel::Coeff);
BOOST_CLASS_EXPORT(bagel::Kinetic);

// Math
BOOST_CLASS_EXPORT(bagel::Matrix);
BOOST_CLASS_EXPORT(bagel::ZMatrix);
BOOST_CLASS_EXPORT(bagel::DIIS<bagel::Matrix>);
BOOST_CLASS_EXPORT(bagel::DIIS<bagel::ZMatrix>);

#endif
