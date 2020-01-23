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

using namespace boost::python;

namespace {

// Converts a std::pair instance to a Python tuple.
template <typename T1, typename T2>
struct std_pair_to_tuple
{
  static PyObject* convert(std::pair<T1, T2> const& p)
  {
    return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr());
  }
  static PyTypeObject const* get_pytype()
  {
    return &PyTuple_Type;
  }
};

// Helper for convenience.
template <typename T1, typename T2>
struct std_pair_to_python_converter
{
  std_pair_to_python_converter()
  {
    boost::python::to_python_converter<std::pair<T1, T2>,
                                       std_pair_to_tuple<T1, T2>,
                                       true // std_pair_to_tuple has get_pytype
                                       >();
  }
};

// Converts a std::vector instance to a Python list.
template <typename T1>
struct std_vector_to_list
{
  static PyObject* convert(std::vector<T1> const& vec)
  {
    boost::python::list list;
    for (auto&& entry : vec) {
      list.append(entry);
    }
    return boost::python::incref(list.ptr());
  }

  static PyTypeObject const* get_pytype()
  {
    return &PyList_Type;
  }
};

template <typename T1>
struct std_vector_to_python_converter
{
  std_vector_to_python_converter()
  {
    boost::python::to_python_converter<std::vector<T1>, std_vector_to_list<T1>, true>();
  }
};


// The iterable_converter is copied from https://stackoverflow.com/a/15940413
// Converts any iterable type from python to C++
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter
{
  /// @note Registers converter from a python iterable type to the
  ///       provided type.
  template <typename Container>
  iterable_converter& from_python()
  {
    boost::python::converter::registry::push_back(&iterable_converter::convertible,
                                                  &iterable_converter::construct<Container>,
                                                  boost::python::type_id<Container>());

    // Support chaining.
    return *this;
  }

  /// @brief Check if PyObject is iterable.
  static void* convertible(PyObject* object)
  {
    return PyObject_GetIter(object) ? object : nullptr;
  }

  /// @brief Convert iterable PyObject to C++ container type.
  ///
  /// Container Concept requirements:
  ///
  ///   * Container::value_type is CopyConstructable.
  ///   * Container can be constructed and populated with two iterators.
  ///     I.e. Container(begin, end)
  template <typename Container>
  static void construct(PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    namespace python = boost::python;
    // Object is a borrowed reference, so create a handle indicting it is
    // borrowed for proper reference counting.
    python::handle<> handle(python::borrowed(object));

    // Obtain a handle to the memory block that the converter has allocated
    // for the C++ type.
    typedef python::converter::rvalue_from_python_storage<Container> storage_type;
    void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

    typedef python::stl_input_iterator<typename Container::value_type> iterator;

    // Allocate the C++ type into the converter's memory block, and assign
    // its handle to the converter's convertible variable.  The C++
    // container is populated by passing the begin and end iterators of
    // the python object to the container's constructor.
    new (storage) Container(iterator(python::object(handle)), // begin
                            iterator()); // end
    data->convertible = storage;
  }
};


} // namespace


template <class Vec>
struct VectorExporter
{
  typedef typename Vec::ScalarType ScalarType;
  typedef typename Dune::XT::LA::VectorInterface<typename Vec::Traits, ScalarType> VectorInterfaceType;
  typedef typename VectorInterfaceType::derived_type derived_type;
  typedef typename Vec::RealType RealType;

  static object buffer(Vec& slf)
  {
    PyObject* py_buf =
        PyMemoryView_FromMemory((char*)&slf[0], slf.size() * sizeof(ScalarType), PyBUF_WRITE); // for python3
    //    PyObject* py_buf = PyBuffer_FromReadWriteMemory(&slf[0], slf.size() * sizeof(ScalarType)); // for python2
    object retval = object(handle<>(py_buf));
    return retval;
  }

  static std::shared_ptr<Vec> create_from_buffer(PyObject* memory_view, const size_t buffer_pos, const size_t vec_size)
  {
    Py_buffer* buffer = PyMemoryView_GET_BUFFER(memory_view);
    ScalarType* cxx_buf = (ScalarType*)buffer->buf;
    return std::make_shared<Vec>(vec_size, cxx_buf + buffer_pos, 0);
  }

  static void export_(const std::string& classname)
  {
    boost::python::type_info info = boost::python::type_id<std::pair<size_t, double>>();
    const boost::python::converter::registration* reg = boost::python::converter::registry::query(info);
    if (reg == nullptr) {
      std_pair_to_python_converter<size_t, double>();
    } else if ((*reg).m_to_python == nullptr) {
      std_pair_to_python_converter<size_t, double>();
    }

    void (Vec::*sub_void)(const derived_type&, derived_type&) const = &Vec::sub;
    derived_type (Vec::*sub_vec)(const derived_type&) const = &Vec::sub;

    void (Vec::*add_void)(const derived_type&, derived_type&) const = &Vec::add;
    derived_type (Vec::*add_vec)(const derived_type&) const = &Vec::add;

    class_<Vec, std::shared_ptr<Vec>>(classname.c_str())
        .def(init<const size_t, const ScalarType, optional<const size_t>>())
        .def("create_from_buffer", &create_from_buffer)
        .staticmethod("create_from_buffer")
        .def("__eq__", &Vec::operator==)
        .def("size", &Vec::size)
        .def("add_to_entry", &Vec::add_to_entry)
        .def("__setitem__", &Vec::set_entry)
        .def("__getitem__", &Vec::get_entry)
        .def("l1_norm", &Vec::l1_norm)
        .def("l2_norm", &Vec::l2_norm)
        .def("sup_norm", &Vec::sup_norm)
        .def("standard_deviation", &Vec::standard_deviation)
        .def("set_all", &Vec::set_all)
        .def("valid", &Vec::valid)
        .add_property("dim", &Vec::size)
        .def("mean", &Vec::mean)
        .def("amax", &Vec::amax)
        .def("sub", sub_void)
        .def("sub", sub_vec)
        .def("add", add_void)
        .def("add", add_vec)
        .def("__add__", add_vec)
        .def("__sub__", sub_vec)
        .def("__iadd__", &Vec::iadd)
        .def("__isub__", &Vec::isub)
        .def("dot", &Vec::dot)
        .def("__mul__", &Vec::dot)
        .def("buffer", &buffer)
        .def("scal", &Vec::scal)
        .def("axpy", &Vec::axpy)
        .def("copy", &Vec::copy);
  }
};


BOOST_PYTHON_MODULE(vectors)
{
  // Vectors
  VectorExporter<Dune::XT::LA::CommonDenseVector<double>>::export_("CommonDenseVector");

  iterable_converter().from_python<std::vector<double>>();
  iterable_converter().from_python<std::vector<size_t>>();

  class_<std::vector<std::vector<Dune::XT::LA::CommonDenseVector<double>>>>("std_vec_of_std_vec_of_la_vecs")
      .def(vector_indexing_suite<std::vector<std::vector<Dune::XT::LA::CommonDenseVector<double>>>>());

  class_<std::vector<Dune::XT::LA::CommonDenseVector<double>>>("std_vec_of_la_vecs")
      .def(vector_indexing_suite<std::vector<Dune::XT::LA::CommonDenseVector<double>>>());

  std_vector_to_python_converter<double>();
  std_vector_to_python_converter<size_t>();
}
