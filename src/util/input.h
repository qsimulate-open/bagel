//insert header

#ifndef __SRC_UTIL_INPUT_H
#define __SRC_UTIL_INPUT_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace bagel {

  // rename as you like
  class PTree {
    protected:
      boost::property_tree::ptree data_;

    public:
      PTree() : data_() {}

      PTree(const boost::property_tree::ptree& i) : data_(i) { }

      PTree(const PTree& o) : data_(o.data_) { }

      PTree(const std::string& input) {
        boost::property_tree::json_parser::read_json(input, data_);
      }

      std::shared_ptr<PTree> get_child(const std::string& key) const {
        return std::make_shared<PTree>(data_.get_child(key));
      }

      std::shared_ptr<PTree> get_child_optional(const std::string& key) const {
        auto out = data_.get_child_optional(key);
        return out ? std::make_shared<PTree>(*out) : std::shared_ptr<PTree>();
      }

      template<typename T> T get(const std::string s) const { return data_.get<T>(s); } 
      template<typename T> T get(const std::string s, const T& t) const { return data_.get<T>(s, t); } 
      template<typename T> void put(const std::string s, const T& o) { data_.put<T>(s, o); } 

      void erase(const std::string key) { data_.erase(key); } 

      // for the time being I use this for debugging 
      // TODO remove
      boost::property_tree::ptree data() const { return data_; }

};

}

#endif
