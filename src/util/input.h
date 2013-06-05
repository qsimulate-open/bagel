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
      PTree(const boost::property_tree::ptree& i) : data_(i) { }

      void read_json(const std::string& input) {
        boost::property_tree::json_parser::read_json(input, data_);
      }

      std::shared_ptr<PTree> get_child(const std::string& key) const {
        return std::make_shared<PTree>(data_.get_child(key));
      }

};

}

#endif
