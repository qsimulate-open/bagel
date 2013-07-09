//
// BAGEL - Parallel electron correlation program.
// Filename: parse.h
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

#ifndef __SRC_INPUT_PARSE_H
#define __SRC_INPUT_PARSE_H

#include <boost/spirit/include/qi.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/bind.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stack>

namespace bagel {

enum class NodeType { base, object, vector };

struct ParseNode {
  NodeType type_;
  std::shared_ptr<boost::property_tree::ptree> data_;

  ParseNode(NodeType type, std::shared_ptr<boost::property_tree::ptree> data) : type_(type), data_(data) {}

  NodeType type() const { return type_; }
  std::shared_ptr<boost::property_tree::ptree> data() const { return data_; }
};

class BagelParser {
  protected:
    std::string filename_;
    std::string contents_;

    std::stack<std::string> key_stack_; // active key stack
    std::stack<ParseNode> node_stack_; // stack of previous nodes up the tree

    std::vector<std::pair<std::string, boost::property_tree::ptree>> base_;

  public:
    BagelParser(std::string filename);

    bool check() const;
    boost::property_tree::ptree parse();

  private:
    void begin_vector();
    void begin_object();

    void close_compound();

    void insert_key(std::string key);
    void insert_value(std::string value);

    template <typename Iterator, typename Skipper>
      struct bagel_checker_grammar : boost::spirit::qi::grammar<Iterator, Skipper> {

        bagel_checker_grammar()
          : bagel_checker_grammar::base_type(bagelinput)
        {
          namespace qi = boost::spirit::qi;

          name = qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9-");

          bagelvalue = +qi::char_("a-zA-Z_0-9+-.");
          bagelvector = qi::lit('[') >> bageldata >> *(qi::lit(',') >> bageldata) >> qi::lit(']');
          bagelobj = qi::lit('{') >> *bagelstatement >> qi::lit('}');

          bageldata = bagelvalue | bagelvector | bagelobj;

          bagelstatement = name >> qi::lit('=') >> bageldata >> qi::lit(';');
          bagelinput = +bagelstatement;
        }

        boost::spirit::qi::rule<Iterator> name;

        boost::spirit::qi::rule<Iterator, Skipper> bagelvalue;
        boost::spirit::qi::rule<Iterator, Skipper> bagelvector;
        boost::spirit::qi::rule<Iterator, Skipper> bagelobj;

        boost::spirit::qi::rule<Iterator, Skipper> bageldata;

        boost::spirit::qi::rule<Iterator, Skipper> bagelstatement;
        boost::spirit::qi::rule<Iterator, Skipper> bagelinput;

    };

    template <typename Iterator, typename Skipper>
      struct bagel_parser_grammar : boost::spirit::qi::grammar<Iterator, Skipper> {

        bagel_parser_grammar(BagelParser* self)
          : bagel_parser_grammar::base_type(bagelinput)
        {
          namespace qi = boost::spirit::qi;

          /* Define Grammar */
          name %= qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9-");

          bagelvalue  %= +qi::char_("a-zA-Z_0-9+-.");

          bagelvector = qi::lit('[') [ boost::bind(&BagelParser::begin_vector, self) ]
                        >> bageldata >> *(qi::lit(',') >> bageldata)
                        >> qi::lit(']') [ boost::bind(&BagelParser::close_compound, self) ];

          bagelobj    = qi::lit('{') [ boost::bind(&BagelParser::begin_object, self) ]
                        >> *bagelstatement
                        >> qi::lit('}') [ boost::bind(&BagelParser::close_compound, self) ];

          bageldata =   bagelvalue [ boost::bind(&BagelParser::insert_value, self, _1) ]
                      | bagelvector
                      | bagelobj;

          bagelstatement = name [ boost::bind(&BagelParser::insert_key, self, _1) ]
                           >> qi::lit('=')
                           >> bageldata
                           >> qi::lit(';');

          bagelinput = +bagelstatement;
        }

        boost::spirit::qi::rule<Iterator, std::string()> name;

        boost::spirit::qi::rule<Iterator, std::string(), Skipper> bagelvalue;
        boost::spirit::qi::rule<Iterator, Skipper> bagelvector;
        boost::spirit::qi::rule<Iterator, Skipper> bagelobj;

        boost::spirit::qi::rule<Iterator, Skipper> bageldata;

        boost::spirit::qi::rule<Iterator, Skipper> bagelstatement;
        boost::spirit::qi::rule<Iterator, Skipper> bagelinput;

    };
};

}

#endif
