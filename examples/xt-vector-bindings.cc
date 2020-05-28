// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <dune/xt/la/container.hh>
#include "xt-vector-bindings.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(vectors)
{
  // Vectors
  using CommonDenseVec = Dune::XT::LA::CommonDenseVector<double>;
  VectorExporter<CommonDenseVec>::export_("CommonDenseVector");

  // from iterable python object to the respective C++ class, allows to call e.g. a C++ function
  // void cpp_function(std::vector<double> vec)
  // with a python list
  iterable_converter().from_python<std::vector<double>>();
  iterable_converter().from_python<std::vector<size_t>>();
  iterable_converter().from_python<std::vector<std::vector<size_t>>>();
  // the other way around, allows to return vectors from C++ functions to python
  std_vector_to_python_converter<double>();
  std_vector_to_python_converter<size_t>();
  std_vector_to_python_converter<std::vector<size_t>>();

  // for some reason, std_vector_to_python_converter<CommonDenseVec> does not work
  class_<std::vector<CommonDenseVec>>("std_vec_of_la_vecs").def(vector_indexing_suite<std::vector<CommonDenseVec>>());

  class_<std::vector<std::vector<CommonDenseVec>>>("std_vec_of_std_vec_of_la_vecs")
      .def(vector_indexing_suite<std::vector<std::vector<CommonDenseVec>>>());
}
