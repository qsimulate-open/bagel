//
// Author : Toru Shiozaki
// Date   : Dec 2011
//

// new interface to the input

#ifndef __NEWINT_UTIL_INPUT_H
#define __NEWINT_UTIL_INPUT_H

#include <map>
#include <list>
#include <string>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

class InputData {
  protected:
    std::list<std::pair<std::string, std::multimap<std::string, std::string> > >data_; 
    const std::string inputfile_;

  public:
    InputData(const std::string filename);
    ~InputData() {};

    std::multimap<std::string, std::string> get_input(const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      if (iter == data_.end())
        throw std::runtime_error(t + " does not appear to be present in your input"); 
      return iter->second;
    };

    bool exist (const std::string t) const {
      auto iter = data_.begin();
      for (; iter != data_.end(); ++iter) if (iter->first == t) break;
      return data_.end() != iter;
    };

    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > data() { return data_; };

};

template <typename T> T read_input(std::multimap<std::string, std::string> idat, const std::string key, const T defvalue) {
  T out = defvalue;
  auto iter = idat.find(key);
  if (iter != idat.end()) out = boost::lexical_cast<T>(iter->second);
  return out;
};

#endif
