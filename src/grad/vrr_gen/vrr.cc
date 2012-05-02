//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include "vrr.h"
#include <fstream>
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
  contents += "#include \"gvrrlist.h\"\n";
  contents += "\n";

  contents += "// returns double array of length " + lexical_cast<string>((a_ + 1) * (c_ + 1) * rank_) + "\n";
  contents += "void GVRRList::" + filename + "(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {\n"; 

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
  string contents;
  int i = 0;
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = 1.0;\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = C00_[" + lexical_cast<string>(t) + "];\n";
  for (int a = 2; a != n + 1; ++a) {
    contents += "\n";
    if (n == 2) {
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_[" + lexical_cast<string>(t) + "];\n";
      }
    } else {
      if (a == 2) contents += "  double B10_current[" + lexical_cast<string>(rank_) + "];\n";
      for (int t = 0; t != rank_; ++t) { 
        contents += "  B10_current[" + lexical_cast<string>(t) + "] " + (a == 2 ? "=" : "+=");
        contents += " B10_[" + lexical_cast<string>(t) + "];\n";
      }
      contents += "\n";
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_current[" + lexical_cast<string>(t) + "]";
        if (a != 2) { 
          contents += " * data_[" + lexical_cast<string>(i - 2 * rank_) + "];\n";
        } else {
          contents += ";\n";
        }
      }
    }
  }
  
  return make_pair(contents, i);
}


const pair<string, int> VRR::vrr0m(const int m) const {
  string contents;
  int i = 0;
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = 1.0;\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = D00_[" + lexical_cast<string>(t) + "];\n";
  for (int c = 2; c != m + 1; ++c) {
    contents += "\n";
    if (m == 2) {
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B01_[" + lexical_cast<string>(t) + "];\n";
      }
    } else {
      if (c == 2) contents += "  double B01_current[" + lexical_cast<string>(rank_) + "];\n";
      for (int t = 0; t != rank_; ++t) { 
        contents += "  B01_current[" + lexical_cast<string>(t) + "] " + (c == 2 ? "=" : "+=");
        contents += " B01_[" + lexical_cast<string>(t) + "];\n";
      }
      contents += "\n";
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B01_current[" + lexical_cast<string>(t) + "]";
        if (c != 2) { 
          contents += " * data_[" + lexical_cast<string>(i - 2 * rank_) + "];\n";
        } else {
          contents += ";\n";
        }
      }
    }
  } 
  return make_pair(contents, i);
}


const pair<string, int> VRR::vrr11() const {
  string contents;
  int i = 0;
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = 1.0;\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = C00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = D00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) {
    contents += "  data_[" + lexical_cast<string>(i) + "]";
    contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
    contents += " + B00_[" + lexical_cast<string>(t) + "];\n"; //* data_[" + lexical_cast<string>(i - rank_ * 3) + "];\n";
  }
  return make_pair(contents, i);
}

/*
const pair<string, int> VRR::vrr1m(const int m) const {
  string contents;
  pair<string, int> zm = vrr0m(m);
  contents += zm.first;
  int i = zm.second;

  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = C00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) { 
    contents += "  data_[" + lexical_cast<string>(i) + "] = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
    contents += " + B00_[" + lexical_cast<string>(t) + "];\n";
  }
  for (int c = 2; c != m + 1; ++c) {
    contents += "\n";
    if (m == 2) {
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B01_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + B00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (m + 2) * rank_) + "];\n";
      }
    } else {
      for (int t = 0; t != rank_; ++t) { 
        contents += "  B01_current[" + lexical_cast<string>(t) + "] " + (c == 2 ? "=" : "+=");
        contents += " B01_[" + lexical_cast<string>(t) + "];\n";
      }
      contents += "\n";
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B01_current[" + lexical_cast<string>(t) + "]";
        contents += " * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + B00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (m + 2) * rank_) + "];\n";
      }
    }
  } 
  return make_pair(contents, i);
}
*/


const pair<string, int> VRR::vrrn1(const int n) const {
  string contents;
  pair<string, int> n0 = vrrn0(n);
  contents = n0.first;
  int i = n0.second;

  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = D00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) { 
    contents += "  data_[" + lexical_cast<string>(i) + "]";
    contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
    contents += " + B00_[" + lexical_cast<string>(t) + "];\n";
  }
  for (int a = 2; a != n + 1; ++a) {
    contents += "\n";
    if (n == 2) {
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + B00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
      }
    } else {
      for (int t = 0; t != rank_; ++t) 
        contents += "  B10_current[" + lexical_cast<string>(t) + "] " + ((a == 2) ? "=" : "+=") + " B10_[" + lexical_cast<string>(t) + "];\n";
      contents += "\n";
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + B00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
      }
    }
  }

  return make_pair(contents, i);
}


const pair<string, int> VRR::vrrnm(const int n, const int m) const {
  string contents;
  pair<string, int> n0 = vrrn0(n);
  contents = n0.first;
  int i = n0.second;

  ////////// c = 1 //////////
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) 
    contents += "  data_[" + lexical_cast<string>(i) + "] = D00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  contents += "  double cB00_current[" + lexical_cast<string>(rank_) + "];\n";
  for (int t = 0; t != rank_; ++t) 
    contents += "  cB00_current[" + lexical_cast<string>(t) + "] = B00_[" + lexical_cast<string>(t) + "];\n";
  contents += "\n";
  for (int t = 0; t != rank_; ++t, ++i) { 
    contents += "  data_[" + lexical_cast<string>(i) + "]";
    contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
    contents += " + cB00_current[" + lexical_cast<string>(t) + "];\n";
  }
  for (int a = 2; a != n + 1; ++a) {
    contents += "\n";
    if (n == 2) {
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + cB00_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
      }
    } else {
      for (int t = 0; t != rank_; ++t) 
        contents += "  B10_current[" + lexical_cast<string>(t) + "] " + ((a == 2) ? "=" : "+=") + " B10_[" + lexical_cast<string>(t) + "];\n";
      contents += "\n";
      for (int t = 0; t != rank_; ++t, ++i) { 
        contents += "  data_[" + lexical_cast<string>(i) + "]";
        contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
        contents += " + B10_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
        contents += " + cB00_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
      }
    }
  }
  ///////////////////////////
  for (int c = 2; c != m + 1; ++c) {
    contents += "\n";
    if (m != 2) {
      if (c == 2)
        contents += "  double B01_current[" + lexical_cast<string>(rank_) + "];\n";
      for (int t = 0; t != rank_; ++t) 
        contents += "  B01_current[" + lexical_cast<string>(t) + "] " + (c == 2 ? "=" : "+=") + " B01_[" + lexical_cast<string>(t) + "];\n";
    }
    contents += "\n";
    for (int t = 0; t != rank_; ++t, ++i) { 
      contents += "  data_[" + lexical_cast<string>(i) + "]";
      contents += " = D00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 1) * rank_) + "]";
      if (m == 2) {
        contents += " + B01_[" + lexical_cast<string>(t) + "];";
      } else {
        contents += " + B01_current[" + lexical_cast<string>(t) + "]";
        contents += c == 2 ? ";" :" * data_[" + lexical_cast<string>(i - 2 * (n + 1) * rank_) + "];";
      } 
      contents += "\n";
    }
    contents += "\n";
    for (int t = 0; t != rank_; ++t) 
      contents += "  cB00_current[" + lexical_cast<string>(t) + "] += B00_[" + lexical_cast<string>(t) + "];\n";
    contents += "\n";
    for (int t = 0; t != rank_; ++t, ++i) { 
      contents += "  data_[" + lexical_cast<string>(i) + "]";
      contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
      contents += " + cB00_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
    }
    for (int a = 2; a != n + 1; ++a) {
      contents += "\n";
      if (n == 2) {
        for (int t = 0; t != rank_; ++t, ++i) { 
          contents += "  data_[" + lexical_cast<string>(i) + "]";
          contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
          contents += " + B10_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
          contents += " + cB00_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
        }
      } else {
        for (int t = 0; t != rank_; ++t) 
          contents += "  B10_current[" + lexical_cast<string>(t) + "] " + ((a == 2) ? "=" : "+=") + " B10_[" + lexical_cast<string>(t) + "];\n";
        contents += "\n";
        for (int t = 0; t != rank_; ++t, ++i) { 
          contents += "  data_[" + lexical_cast<string>(i) + "]";
          contents += " = C00_[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - rank_) + "]";
          contents += " + B10_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - 2 * rank_) + "]";
          contents += " + cB00_current[" + lexical_cast<string>(t) + "] * data_[" + lexical_cast<string>(i - (n + 2) * rank_) + "];\n";
        }
      }
    }
  }

  return make_pair(contents, i);
}


