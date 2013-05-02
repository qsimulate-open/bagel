//
// BAGEL - Parallel electron correlation program.
// Filename: test_main.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites

#include <stddef.h>
#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <boost/test/unit_test.hpp>
#include <src/parallel/resources.h>
#include <src/parallel/mpi_interface.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace bagel;

Resources b(8);
Resources* bagel::resources__ = &b;
static int argc = 0;
static char** argv;
MPI_Interface c(argc, argv);
MPI_Interface* bagel::mpi__ = &c;

static double THRESH = 1.0e-8;

bool compare(const double a, const double b, const double thr = THRESH) { return fabs(a-b) < thr; };

template<class T>
bool compare(const T a, const T b, const double thr = THRESH) {
 if (a.size() != b.size()) throw std::logic_error("comparing vectors with different sizes");
 bool out = true;
 for (auto i = a.begin(), j = b.begin(); i != a.end(); ++i, ++j) out &= compare(*i, *j, thr);
 return out;
};

#include <src/scf/test_scf.cc>

#include <src/ks/test_ks.cc>

#include <src/rel/test_rel.cc>

#include <src/prop/test_prop.cc>

#include <src/mp2/test_mp2.cc>

#include <src/casscf/test_casscf.cc>

#include <src/fci/test_fci.cc>

#include <src/opt/test_opt.cc>

#if 0
#include <src/io/test_molden.cc>

#include <src/util/test_localize.cc>
#endif

#include <src/meh/test_meh.cc>
