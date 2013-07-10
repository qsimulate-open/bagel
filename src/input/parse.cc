//
// BAGEL - Parallel electron correlation program.
// Filename: parse.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/input/parse.h>

using namespace std;
using namespace bagel;

using ptree = boost::property_tree::ptree;

BagelParser::BagelParser(string filename) : filename_(filename) {
  ifstream in(filename_, ios_base::in);

  if (!in.good()) throw runtime_error(string("Error! Could not open input file: ") + filename_);

  in.unsetf(std::ios::skipws); // No white space skipping!
  contents_ = string( (std::istream_iterator<char>(in)) , (std::istream_iterator<char>()) );
}

bool BagelParser::check() const {
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;

  // Set iterator type
  using Iterator = string::const_iterator;

  Iterator iter = contents_.begin();
  Iterator end = contents_.end();

  // Skip white-space and comments
  qi::rule<Iterator> skipper = ( ascii::space
                         | (qi::lit("/*") >> *(qi::char_ - qi::lit("*/")) >> qi::lit("*/"))
                         | (qi::lit("//") >> *(qi::char_ - qi::eol) >> qi::eol)
                         );

  // Instantiate checker grammar
  bagel_checker_grammar<Iterator, decltype(skipper)> checker;

  // Parse
  bool result = qi::phrase_parse(iter, end, checker, skipper);

  return (iter == end);
}

ptree BagelParser::parse() {
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;

  // Set iterator type
  using Iterator = string::const_iterator;

  // Setup necessary parts of the parser
  Iterator iter = contents_.begin();
  Iterator end = contents_.end();

  // Skip white-space and comments
  qi::rule<Iterator> skipper = ( ascii::space
                         | (qi::lit("/*") >> *(qi::char_ - qi::lit("*/")) >> qi::lit("*/"))
                         | (qi::lit("//") >> *(qi::char_ - qi::eol) >> qi::eol)
                         );

  // Instantiate parser grammar
  bagel_parser_grammar<Iterator, decltype(skipper)> parser(this);

  base_.clear();
  while (!node_stack_.empty()) node_stack_.pop();
  node_stack_.push(ParseNode(NodeType::base, nullptr));

  // parse
  bool result;
  try {
    result = qi::phrase_parse(iter, end, parser, skipper);
  }
  catch (const qi::expectation_failure<Iterator>& e) {
    cerr << "Error parsing file \'" << filename_ << "\'. Expected " << e.what_.tag << " " << get_error_position(contents_.cbegin(), e.first - 1) << endl;
    throw runtime_error("Failed parsing input file!");
  }
  if ( (!result) || (iter != end) ) throw bagel_parser_error("Failed parsing input file, probably incorrect format.");

  // collect into the appropriate ptree
  ptree bageltree;
  for (auto& ibase : base_) {
    ibase.second.put("title", ibase.first);
    bageltree.push_back(make_pair("",ibase.second));
  }

  ptree out;
  out.add_child("bagel", bageltree);

  return out;
}

void BagelParser::begin_vector() { node_stack_.push(ParseNode(NodeType::vector, make_shared<ptree>())); }

void BagelParser::begin_object() { node_stack_.push(ParseNode(NodeType::object, make_shared<ptree>())); }

void BagelParser::close_compound() {
  auto finished = node_stack_.top();
  node_stack_.pop();

  auto current = node_stack_.top();

  if (current.type() == NodeType::vector) {
    current.data()->push_back(make_pair("", *finished.data()));
  }
  else if (current.type() == NodeType::object) {
    current.data()->add_child(key_stack_.top(), *finished.data());
    key_stack_.pop();
  }
  else /* type() == NodeType::base */ {
    base_.push_back(make_pair(key_stack_.top(), *finished.data()));
    key_stack_.pop();
  }
}

void BagelParser::insert_value(string value) {
  auto current = node_stack_.top();
  if (current.type() == NodeType::object) {
    current.data()->put(key_stack_.top(), value);
    key_stack_.pop();
  }
  else if (current.type() == NodeType::vector) {
    ptree child;
    child.put("", value);
    current.data()->push_back(make_pair("", child));
  }
  else /* type() == NodeType::base */ {
    ptree child;
    child.put_value(value);
    base_.push_back(make_pair(key_stack_.top(), child));
    key_stack_.pop();
  }
}

void BagelParser::insert_key(string key) { key_stack_.push(key); }

string BagelParser::get_error_position(std::string line) {
  int line_number = 1;
  std::size_t found = line.find('\n');
  while(found != std::string::npos) {
    ++line_number;
    found = line.find('\n', found+1);
  }

  int col_number;
  std::size_t last = line.find_last_of('\n');
  if(last == std::string::npos) col_number = line.size();
  else col_number = line.size() - (last + 1);

  std::stringstream ss;
  ss << "on line " << line_number << ", at position " << col_number << ".";
  return ss.str();
}
