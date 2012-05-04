//
// Newint - Parallel electron correlation program.
// Filename: vrr.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include "vrr.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

const string VRR::dump(const string filename) const {
  string header;
  string contents;

  contents += "//\n";
  contents += "// Author: Toru Shiozaki\n";
  contents += "// Machine Generated Code in NewInt\n";
  contents += "//\n";
  contents += "\n";
  contents += "#include \"vrrlist.h\"\n";
  contents += "\n";

  contents += "// returns double array of length " + lexical_cast<string>((a_ + 1) * (c_ + 1) * rank_) + "\n";
  contents += "void VRRList::" + filename + "(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {\n"; 

  if (a_ > 1 && c_ > 1)        contents += vrrnm(a_, c_).first; 
  else if (a_ == 1 && c_ >  1) contents += vrrnm(1 , c_).first;
  else if (a_ == 0 && c_ >  0) contents += vrr0m(c_).first;
  else if (a_ >  1 && c_ == 1) contents += vrrn1(a_).first;
  else if (a_ == 1 && c_ == 1) contents += vrr11(  ).first;
  else if (a_ >  0 && c_ == 0) contents += vrrn0(a_).first; 
  else if (a_ == 0 && c_ == 0) contents += vrr00(  ).first;
  else cout << "NEVER COMES HERE" << endl;
    

  contents += "}\n";

  ofstream cc;
  const string ccfilename = filename + ".cc";
  cc.open(ccfilename.c_str());
  cc << contents << endl;
  cc.close();

  return header; 
}


const pair<string, int> VRR::vrr00() const {
  string contents;
  int i = 0;
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = 1.0;\n";
  return make_pair(contents, i);
}


const pair<string, int> VRR::vrrn0(const int n) const {
  stringstream contents;
  int i = 0;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = 1.0;\n";
  i += rank_;
  contents << endl;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = C00_[t];" << endl;
  i += rank_;
  for (int a = 2; a != n + 1; ++a) {
    contents << endl;
    if (n == 2) {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B10_[t];\n";
      i += rank_;
    } else {
      if (a == 2) contents << "  double B10_current[" << rank_ << "];\n";
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
      contents << "    B10_current[t] " << (a == 2 ? "=" : "+=");
      contents << " B10_[t];" << endl;
      contents << endl;
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B10_current[t]";
      if (a != 2) { 
        contents << " * data_[" << (i - 2 * rank_) << "+t];";
      } else {
        contents << ";";
      }
      contents << endl;
      i += rank_;
    }
  }
  
  return make_pair(contents.str(), i);
}


const pair<string, int> VRR::vrr0m(const int m) const {
  stringstream contents;
  int i = 0;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = 1.0;\n";
  contents << endl;
  i += rank_;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = D00_[t];\n";
  i += rank_;
  for (int c = 2; c != m + 1; ++c) {
    contents << endl;
    if (m == 2) {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = D00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B01_[t];" << endl;
      i += rank_;
    } else {
      if (c == 2) contents << "  double B01_current[" << rank_ << "];\n";
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    B01_current[t] " << (c == 2 ? "=" : "+=");
      contents << " B01_[t];" << endl;
      contents << endl;
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = D00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B01_current[t]";
      if (c != 2) { 
        contents << " * data_[" << (i - 2 * rank_) << "+t];";
      } else {
        contents << ";";
      }
      contents << endl;
      i += rank_;
    }
  } 
  return make_pair(contents.str(), i);
}


const pair<string, int> VRR::vrr11() const {
  stringstream contents;
  int i = 0;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = 1.0;\n";
  i += rank_;
  contents << "\n";
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = C00_[t];\n";
  contents << "\n";
  i += rank_;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
  contents << "    data_[" << i << "+t] = D00_[t];\n";
  contents << "\n";
  i += rank_;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
  contents << "    data_[" << i << "+t]";
  contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
  contents << " + B00_[t];\n"; //* data_[" + lexical_cast<string>(i - rank_ * 3) + "];\n";
  return make_pair(contents.str(), i);
}


const pair<string, int> VRR::vrrn1(const int n) const {
  stringstream contents;
  pair<string, int> n0 = vrrn0(n);
  contents << n0.first;
  int i = n0.second;

  contents << "\n";
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = D00_[t];\n";
  contents << "\n";
  i += rank_;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
  contents << "    data_[" << i << "+t]";
  contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
  contents << " + B00_[t];\n";
  i += rank_;
  for (int a = 2; a != n + 1; ++a) {
    contents << "\n";
    if (n == 2) {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B10_[t] * data_[" << (i - 2 * rank_) << "+t]";
      contents << " + B00_[t] * data_[" << (i - (n + 2) * rank_) << "+t];\n";
      i += rank_;
    } else {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    B10_current[t] " << ((a == 2) ? "=" : "+=") << " B10_[t];\n";
      contents << endl;
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B10_current[t] * data_[" << (i - 2 * rank_) << "+t]";
      contents << " + B00_[t] * data_[" << (i - (n + 2) * rank_) << "+t];\n";
      i += rank_;
    }
  }

  return make_pair(contents.str(), i);
}


const pair<string, int> VRR::vrrnm(const int n, const int m) const {
  stringstream contents;
  pair<string, int> n0 = vrrn0(n);
  contents << n0.first;
  int i = n0.second;

  ////////// c = 1 //////////
  contents << endl;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
  contents << "    data_[" << i << "+t] = D00_[t];" << endl;
  i += rank_;
  contents << endl;
  contents << "  double cB00_current[" << rank_ << "];\n";
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
  contents << "    cB00_current[t] = B00_[t];" << endl;
  contents << endl;
  contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
  contents << "    data_[" << i << "+t] = C00_[t] * data_[" << i - rank_ << "+t] + cB00_current[t];" << endl;
  i += rank_;

  for (int a = 2; a != n + 1; ++a) {
    contents << endl;
    if (n == 2) {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << i - rank_ << "+t]";
      contents << " + B10_[t] * data_[" << i - 2 * rank_ << "+t]";
      contents << " + cB00_current[t] * data_[" << (i - (n + 2) * rank_) << "+t];\n";
      i += rank_;
    } else {
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    B10_current[t] " << ((a == 2) ? "=" : "+=") << " B10_[t];\n";
      contents << endl;
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    data_[" << i << "+t]";
      contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
      contents << " + B10_current[t] * data_[" << (i - 2 * rank_) << "+t]";
      contents << " + cB00_current[t] * data_[" << (i - (n + 2) * rank_) << "+t];\n";
      i += rank_;
    }
  }
  ///////////////////////////
  for (int c = 2; c != m + 1; ++c) {
    contents << "\n";
    if (m != 2) {
      if (c == 2)
      contents << "  double B01_current[" << rank_ << "];" << endl;
      contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
      contents << "    B01_current[t] " << (c == 2 ? "=" : "+=") << " B01_[t];\n";
    }
    contents << endl;
    contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
    contents << "    data_[" << i << "+t]";
    contents << " = D00_[t] * data_[" << (i - (n + 1) * rank_) << "+t]";
    if (m == 2) {
      contents << " + B01_[t];";
    } else {
      contents << " + B01_current[t]" << (c == 2 ? ";" :" * data_[" + lexical_cast<string>(i-2*(n+1)*rank_) + "+t];");
    } 
    contents << endl;
    i += rank_;
    contents << endl;
    contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
    contents << "    cB00_current[t] += B00_[t];\n";
    contents << endl;
    contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
    contents << "    data_[" << i << "+t]";
    contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
    contents << " + cB00_current[t] * data_[" << (i - (n + 2) * rank_) << "+t];" << endl;
    i += rank_;
    for (int a = 2; a != n + 1; ++a) {
      contents << endl;
      if (n == 2) {
        contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
        contents << "    data_[" << i << "+t]";
        contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
        contents << " + B10_[t] * data_[" << (i - 2 * rank_) << "+t]";
        contents << " + cB00_current[t] * data_[" << (i - (n + 2) * rank_) << "+t];" << endl;
        i += rank_;
      } else {
        contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl; 
        contents << "    B10_current[t] " << ((a == 2) ? "=" : "+=") << " B10_[t];" << endl;
        contents << endl;
        contents << "  for (int t = 0; t != " << rank_ << "; ++t)" << endl;
        contents << "    data_[" << i << "+t]";
        contents << " = C00_[t] * data_[" << (i - rank_) << "+t]";
        contents << " + B10_current[t] * data_[" << (i - 2 * rank_) << "+t]";
        contents << " + cB00_current[t] * data_[" << (i - (n + 2) * rank_) << "+t];" << endl;
        i += rank_;
      }
    }
  }

  return make_pair(contents.str(), i);
}


