//copy header

#include <src/util/input.h>
#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(doublewrap){
  class_<bagel::PTree>("PTree", init<boost::property_tree::ptree>())
  .def("get_child", &bagel::PTree::get_child)
;
}
