//
// BAGEL - Parallel electron correlation program.
// Filename: input.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <src/util/input.h>
#include <boost/regex.hpp>

using namespace std;
using namespace bagel;

static vector<string> split(const string& s, const string& delimiters) {
  vector<string> out;
  size_t current;
  size_t next = -1;
  do {
    current = next + 1;
    next = s.find_first_of(delimiters, current);
    const string now = s.substr(current, next - current);
    out.push_back(now);
  } while (next != std::string::npos);
  return out;
}
static string &ltrim(string &s) { s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace)))); return s; }
static string &rtrim(string &s) { s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end()); return s; }
static string &trim(string &s) { return ltrim(rtrim(s)); }

InputData::InputData(const string filename) : inputfile_(filename) {
  ifstream ifs;
  ifs.open(inputfile_);
  if (!ifs.is_open())
    throw runtime_error("Your input file cannot be opened.");

  string content = "";
  while(1) {
    string sline;
    if (!getline(ifs, sline)) break;
    // get rid of comments
    const boost::regex comm("^(.*)//(.*)$");
    const string ss = sline;
    auto start = ss.begin();
    auto end = ss.end();
    boost::smatch what;
    while (boost::regex_search(start, end, what, comm)) {
      sline = what[1];
      const string sk = sline;
      start = sk.begin();
      end   = sk.end();
    }
    content += sline;
  }
  // true false to 1 and 0 to be friendly to lexical_cast which does not accept true/false
  {
    const boost::regex reg_true("^(.*)=(\\s*)(true)(\\s*);(.*)$");
    const boost::regex reg_false("^(.*)=(\\s*)(false)(\\s*);(.*)$");
    const string ss = content;
    auto start = ss.begin();
    auto end = ss.end();
    boost::smatch what;
    while (boost::regex_search(start, end, what, reg_true)) {
      stringstream ss; ss << what[1] << "=1;" << what[5];
      const string sk = ss.str();
      start = sk.begin();
      end = sk.end();
      content = sk;
    }
    while (boost::regex_search(start, end, what, reg_false)) {
      stringstream ss; ss << what[1] << "=0;" << what[5];
      const string sk = ss.str();
      start = sk.begin();
      end = sk.end();
      content = sk;
    }
  }

  // first split with { and }
  if (content.find("}") == string::npos) throw runtime_error("Check your input format");
  vector<string> blocks = split(content, "}");
  for (auto& it : blocks) {
    if (it.size() == 0) continue;
    vector<string> pairs = split(it, "{");

    trim(pairs[0]);
    string tag = pairs[0];
    transform(tag.begin(), tag.end(), tag.begin(), ::tolower);

    vector<string> lines = split(pairs[1], ";");
    multimap<string, string> tmp;
    for (auto& i : lines) {
      if (i.size() > 0) {
        vector<string> ll = split(i, "=");
        if (ll.size() != 2) {
          stringstream ss;
          ss << "There seem " << ll.size()-1 << " '=' in a single directive. Check your input.";
          throw runtime_error(ss.str());
        }
        trim(ll[0]); trim(ll[1]);
        transform(ll[0].begin(), ll[0].end(), ll[0].begin(), ::tolower);
        transform(ll[1].begin(), ll[1].end(), ll[1].begin(), ::tolower);
        tmp.insert(make_pair(ll[0],ll[1]));
      }
    }
    data_.push_back(make_pair(tag, tmp));
  }
}
