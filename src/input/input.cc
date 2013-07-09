//
// BAGEL - Parallel electron correlation program.
// Filename: input.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/input/input.h>
#include <src/input/parse.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
//#include <iostream>

using namespace bagel;
using namespace std;

PTree::PTree(const std::string& input) {
  // Check first from clues from file extension
  if (string(input.end() - 5, input.end()) == ".json") {
    boost::property_tree::json_parser::read_json(input, data_);
  }
  else if (string(input.end() - 4, input.end()) == ".xml") {
    boost::property_tree::xml_parser::read_xml(input, data_);
  }
  else if (string(input.end() - 4, input.end()) == ".bgl" || string(input.end() - 6, input.end()) == ".bagel"){
    BagelParser bp(input);
    data_ = bp.parse();
  }
  else { // Unhelpful file extension -> just try them all!
    try {
      boost::property_tree::json_parser::read_json(input, data_);
    }
    catch (boost::property_tree::json_parser_error& e) {
      try {
        boost::property_tree::xml_parser::read_xml(input, data_);
      }
      catch (boost::property_tree::xml_parser_error& f) {
        try {
          BagelParser bp(input);
          data_ = bp.parse();
        }
        catch (bagel_parser_error& g) {
          throw runtime_error("Failed to determine input file format. Try specifying it with the file extension ( \'.json\', \'.xml\', or (\'.bgl\'|\'.bagel\' ).");
        }
      }
    }

  }
}


//Dereference operator - return the current node's data.
const std::shared_ptr<const PTree> PTreeIterator::operator*() { return make_shared<const PTree>(current->second); }

PTreeIterator PTree::begin() const { return PTreeIterator(data_.begin()); }
PTreeIterator PTree::end()   const { return PTreeIterator(data_.end());   }


size_t PTree::size() const {
  return data_.size();
}


void PTree::print() const {
  write_json(cout, data_);
}
