//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: input.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <fstream>
#include <string>
#include <src/util/input/input.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::PTree)

using namespace bagel;
using namespace std;

PTree::PTree(const string& input) {
  // Check first from clues from file extension
  const size_t n = input.find_last_of(".");
  const string extension = (n != string::npos) ? input.substr(n) : "";
  if (extension == ".json") {
    try {
      boost::property_tree::json_parser::read_json(input, data_);

    // making error messages more user-friendly
    } catch (const exception& e) {
      string errstring(e.what());
      size_t p1 = errstring.find("(");
      size_t p2 = errstring.find(")");
      string error_message = "Syntax error near line " + errstring.substr(p1+1, p2-p1-1) + " of " + errstring.substr(0, p1) + ".  ";
      if (errstring.find("expected ']' or ','") != string::npos || errstring.find("expected '}' or ','") != string::npos)
        throw runtime_error(error_message + "Probably you forgot a comma or bracket.");
      else if (errstring.find("expected value") != string::npos)
        throw runtime_error(error_message + "Expected a value but did not find one...  probably you have an extra comma.");
      else if (errstring.find("expected key string") != string::npos)
        throw runtime_error(error_message + "Expected a keyword but did not find one...  probably you have an extra comma.");
      else
        throw;
    }
  } else if (extension == ".xml") {
    boost::property_tree::xml_parser::read_xml(input, data_);
  } else { // Unhelpful file extension -> just try them all!
    try {
      boost::property_tree::json_parser::read_json(input, data_);
    } catch (boost::property_tree::json_parser_error& e) {
      try {
        boost::property_tree::xml_parser::read_xml(input, data_);
      } catch (boost::property_tree::xml_parser_error& g) {
        throw runtime_error("Failed to determine input file format. Try specifying it with the file extension ( \'.json\' or \'.xml\').");
      }
    }
  }
}


PTree::PTree(stringstream& input) {
  try {
    boost::property_tree::json_parser::read_json(input, data_);
  } catch (boost::property_tree::json_parser_error& e) {
    throw runtime_error("Failed to recognize the input string as a stringify'ed json object");
  }
}


//Dereference operator - return the current node's data.
const shared_ptr<const PTree> PTreeIterator::operator*() { return make_shared<const PTree>(current->second, current->first); }
const shared_ptr<const PTree> PTreeReverseIterator::operator*() { return make_shared<const PTree>(current->second, current->first); }

PTreeIterator PTree::begin() const { return PTreeIterator(data_.begin()); }
PTreeIterator PTree::end()   const { return PTreeIterator(data_.end());   }
PTreeReverseIterator PTree::rbegin() const { return PTreeReverseIterator(data_.rbegin()); }
PTreeReverseIterator PTree::rend()   const { return PTreeReverseIterator(data_.rend());   }

namespace bagel {
template <> void PTree::push_back<shared_ptr<PTree>>(const shared_ptr<PTree>& pt) {
  data_.push_back({"", pt->data_});
}
}


size_t PTree::size() const {
  return data_.size();
}


void PTree::print() const {
  write_json(cout, data_);
}


shared_ptr<const PTree> PTree::read_basis(string name) {
  // convert name to lowercase so things like cc-pVDZ are read
  const int split = name.find_last_of("/");
  name = name.substr(0, split+1) + to_lower(name.substr(split+1));

  shared_ptr<const PTree> out;
  // first try the absolute path (or current directory)
  try {
    out = make_shared<const PTree>(name);
  } catch (...) {
    // next, the standard install location
    try {
      const string prefix(BASIS_DIR);
      const string filename = prefix + "/" + name + ".json";
      out = make_shared<const PTree>(filename);
    } catch (...) {
      // last, the debug location
      const string filename = "../../src/basis/" + name + ".json";
      try {
        out = make_shared<const PTree>(filename);
      } catch (...) {
        throw runtime_error(name + " cannot be opened. Please see if the file is in " + string(BASIS_DIR) + ".\n "
                                 + " You can also specify the full path to the basis file.");
      }
    }
  }
  return out;
}


shared_ptr<PTree> PTree::get_child(const string& key) const {
  auto out = data_.get_child_optional(key);
  if (!out) {
    string thiskey = get<string>("title", key_);
    throw runtime_error("A required keyword is missing from the input" + (thiskey.size() ? " for " + thiskey : "") + ":  " + key);
  }
  return make_shared<PTree>(*out, key);
}


shared_ptr<PTree> PTree::get_child_optional(const string& key) const {
  auto out = data_.get_child_optional(key);
  return out ? make_shared<PTree>(*out, key) : nullptr;
}
