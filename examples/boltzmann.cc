// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

#include "cellmodel.hh"
#include "boltzmann.hh"

using BoltzmannSolver2d = BoltzmannSolver<2>;
using BoltzmannSolver3d = BoltzmannSolver<3>;

// Python bindings
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

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
  typedef typename XT::LA::VectorInterface<typename Vec::Traits, ScalarType> VectorInterfaceType;
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


#include <dune/xt/common/disable_warnings.hh>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads2d, BoltzmannSolver2d::init, 0, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_timesteps_overloads2d, BoltzmannSolver2d::next_n_timesteps, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads2d, BoltzmannSolver2d::apply_rhs_operator, 3, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(init_overloads3d, BoltzmannSolver3d::init, 0, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(next_n_timesteps_overloads3d, BoltzmannSolver3d::next_n_timesteps, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(apply_rhs_overloads3d, BoltzmannSolver3d::apply_rhs_operator, 3, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(visualize_overloads, CellModelSolver::visualize, 3, 6)
#include <dune/xt/common/reenable_warnings.hh>


BOOST_PYTHON_MODULE(libhapodgdt)
{
  using VectorType = typename BoltzmannSolver3d::VectorType;
  using RangeFieldType = typename BoltzmannSolver3d::RangeFieldType;

  // 2d
  VectorType (BoltzmannSolver2d::*apply_rhs_without_params2d)(VectorType, const double) const =
      &BoltzmannSolver2d::apply_rhs_operator;
  VectorType (BoltzmannSolver2d::*apply_rhs_with_params2d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver2d::apply_rhs_operator;

  class_<BoltzmannSolver2d>("BoltzmannSolver2d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double,
                                          const bool>>())
      .def("init", &BoltzmannSolver2d::init, init_overloads2d())
      .def("linear", &BoltzmannSolver2d::linear)
      .def("solve", &BoltzmannSolver2d::solve)
      .def("next_n_timesteps", &BoltzmannSolver2d::next_n_timesteps, next_n_timesteps_overloads2d())
      .def("reset", &BoltzmannSolver2d::reset)
      .def("finished", &BoltzmannSolver2d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver2d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver2d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver2d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver2d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params2d)
      .def("apply_rhs_operator", apply_rhs_with_params2d, apply_rhs_overloads2d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver2d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver2d::get_initial_values)
      .def("current_time", &BoltzmannSolver2d::current_time)
      .def("set_current_time", &BoltzmannSolver2d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver2d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver2d::time_step_length)
      .def("t_end", &BoltzmannSolver2d::t_end);

  // 3d
  VectorType (BoltzmannSolver3d::*apply_rhs_without_params3d)(VectorType, const double) const =
      &BoltzmannSolver3d::apply_rhs_operator;
  VectorType (BoltzmannSolver3d::*apply_rhs_with_params3d)(VectorType,
                                                           const double,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType,
                                                           const RangeFieldType) =
      &BoltzmannSolver3d::apply_rhs_operator;

  class_<BoltzmannSolver3d>("BoltzmannSolver3d",
                            init<optional<const std::string,
                                          const size_t,
                                          const size_t,
                                          const bool,
                                          const bool,
                                          const double,
                                          const double,
                                          const double,
                                          const double,
                                          const bool>>())
      .def("init", &BoltzmannSolver3d::init, init_overloads3d())
      .def("linear", &BoltzmannSolver3d::linear)
      .def("solve", &BoltzmannSolver3d::solve)
      .def("next_n_time_steps", &BoltzmannSolver3d::next_n_timesteps, next_n_timesteps_overloads3d())
      .def("reset", &BoltzmannSolver3d::reset)
      .def("finished", &BoltzmannSolver3d::finished)
      .def("apply_kinetic_operator", &BoltzmannSolver3d::apply_kinetic_operator)
      .def("apply_restricted_kinetic_operator", &BoltzmannSolver3d::apply_restricted_kinetic_operator)
      .def("prepare_restricted_operator", &BoltzmannSolver3d::prepare_restricted_operator)
      .def("source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs)
      .def("len_source_dofs", &BoltzmannSolver3d::restricted_op_input_dofs_size)
      .def("apply_rhs_operator", apply_rhs_without_params3d)
      .def("apply_rhs_operator", apply_rhs_with_params3d, apply_rhs_overloads3d())
      .def("set_rhs_operator_parameters", &BoltzmannSolver3d::set_rhs_operator_parameters)
      .def("get_initial_values", &BoltzmannSolver3d::get_initial_values)
      .def("current_time", &BoltzmannSolver3d::current_time)
      .def("set_current_time", &BoltzmannSolver3d::set_current_time)
      .def("set_current_solution", &BoltzmannSolver3d::set_current_solution)
      .def("time_step_length", &BoltzmannSolver3d::time_step_length)
      .def("t_end", &BoltzmannSolver3d::t_end);

  // Cell Model
  // Boost.Python cannot handle more than about 14 arguments in the constructor by default
  // TODO: Find a way to increase that limit in boost or change constructor.
  class_<CellModelSolver>("CellModelSolver",
                          init<optional<const std::string,
                                        const double,
                                        const unsigned int,
                                        const unsigned int,
                                        const bool,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        const double,
                                        // const double,
                                        // const double,
                                        // const double,
                                        const double>>())
      .def("visualize", &CellModelSolver::visualize, visualize_overloads())
      .def("prepare_pfield_operator", &CellModelSolver::prepare_pfield_operator)
      .def("prepare_ofield_operator", &CellModelSolver::prepare_ofield_operator)
      .def("prepare_stokes_operator", &CellModelSolver::prepare_stokes_operator)
      .def("apply_pfield_operator", &CellModelSolver::apply_pfield_operator)
      .def("apply_ofield_operator", &CellModelSolver::apply_ofield_operator)
      .def("apply_stokes_operator", &CellModelSolver::apply_stokes_operator)
      .def("set_pfield_variables", &CellModelSolver::set_pfield_variables)
      .def("set_ofield_variables", &CellModelSolver::set_ofield_variables)
      .def("set_stokes_variables", &CellModelSolver::set_stokes_variables)
      .def("solve", &CellModelSolver::solve)
      .def("next_n_timesteps", &CellModelSolver::next_n_timesteps)
      .def("solve_pfield", &CellModelSolver::solve_pfield)
      .def("solve_ofield", &CellModelSolver::solve_ofield)
      .def("solve_stokes", &CellModelSolver::solve_stokes)
      .def("apply_pfield_product_operator", &CellModelSolver::apply_pfield_product_operator)
      .def("apply_ofield_product_operator", &CellModelSolver::apply_ofield_product_operator)
      .def("apply_stokes_product_operator", &CellModelSolver::apply_stokes_product_operator)
      .def("num_cells", &CellModelSolver::num_cells)
      .def("finished", &CellModelSolver::finished)
      .def("linear", &CellModelSolver::linear)
      .def("pfield_vector", &CellModelSolver::pfield_vector)
      .def("ofield_vector", &CellModelSolver::ofield_vector)
      .def("stokes_vector", &CellModelSolver::stokes_vector);

  // Vectors
  VectorExporter<typename XT::LA::CommonDenseVector<double>>::export_("CommonDenseVector");

  iterable_converter().from_python<std::vector<double>>();
  iterable_converter().from_python<std::vector<size_t>>();

  class_<std::vector<std::vector<typename XT::LA::CommonDenseVector<double>>>>("std_vec_of_std_vec_of_la_vecs")
      .def(vector_indexing_suite<std::vector<std::vector<typename XT::LA::CommonDenseVector<double>>>>());

  class_<std::vector<typename XT::LA::CommonDenseVector<double>>>("std_vec_of_la_vecs")
      .def(vector_indexing_suite<std::vector<typename XT::LA::CommonDenseVector<double>>>());

  std_vector_to_python_converter<double>();
  std_vector_to_python_converter<size_t>();
}


int main(int argc, char* argv[])
{
  static const size_t dimDomain = 2;
  static const bool linear = true;
  try {
    // parse options
    if (argc == 1)
      std::cout << "The following options are available: " << argv[0]
                << " [-output_dir DIR -num_save_steps INT -gridsize INT "
                << "  -sigma_s_1 FLOAT -sigma_s_2 FLOAT -sigma_a_1 FLOAT -sigma_a_2 FLOAT"
                << " --no_visualization --silent --random_parameters]" << std::endl;

    size_t num_save_steps = 10;
    size_t grid_size = 20;
    bool visualize = true;
    bool silent = false;
    bool random_parameters = false;
    bool parameters_given = false;
    std::string output_dir;
    double sigma_s_lower = 0, sigma_s_upper = 8, sigma_a_lower = 0, sigma_a_upper = 8;
    double sigma_s_1 = 1, sigma_s_2 = 0, sigma_a_1 = 0, sigma_a_2 = 10;
    for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-output_dir") {
        if (i + 1 < argc) {
          output_dir = argv[++i];
        } else {
          std::cerr << "-output_dir option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-num_save_steps") {
        if (i + 1 < argc) {
          num_save_steps = XT::Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-num_save_steps option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "--no_visualization") {
        visualize = false;
      } else if (std::string(argv[i]) == "--silent") {
        silent = true;
      } else if (std::string(argv[i]) == "--random_parameters") {
        if (parameters_given) {
          std::cerr << "You specified a value for at least one parameter so you can't use --random_parameters!"
                    << std::endl;
          return 1;
        }
        random_parameters = true;
        RandomNumberGeneratorType rng{std::random_device()()};
        std::uniform_real_distribution<double> sigma_s_dist(sigma_s_lower, sigma_s_upper);
        std::uniform_real_distribution<double> sigma_a_dist(sigma_a_lower, sigma_a_upper);
        sigma_s_1 = sigma_s_dist(rng);
        sigma_s_2 = sigma_s_dist(rng);
        sigma_a_1 = sigma_a_dist(rng);
        sigma_a_2 = sigma_a_dist(rng);
      } else if (std::string(argv[i]) == "-gridsize") {
        if (i + 1 < argc) {
          grid_size = XT::Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-gridsize option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_1 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_2 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_2 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_1 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_2 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_2 option requires one argument." << std::endl;
          return 1;
        }
      } else {
        std::cerr << "Unknown option " << std::string(argv[i]) << std::endl;
        return 1;
      }
    }

    std::ofstream parameterfile;
    parameterfile.open(output_dir + "_parameters.txt");
    parameterfile << "Gridsize: " << XT::Common::to_string(grid_size) + " x " + XT::Common::to_string(grid_size)
                  << std::endl;

    // run solver
    parameterfile << "Domain was composed of two materials, parameters were: " << std::endl
                  << "First material: sigma_s = " + XT::Common::to_string(sigma_s_1)
                         + ", sigma_a = " + XT::Common::to_string(sigma_a_1)
                  << std::endl
                  << "Second material: sigma_s = " + XT::Common::to_string(sigma_s_2)
                         + ", sigma_a = " + XT::Common::to_string(sigma_a_2)
                  << std::endl;

    auto solver = std::make_shared<BoltzmannSolver<dimDomain>>(
        output_dir, num_save_steps, grid_size, visualize, silent, sigma_s_1, sigma_s_2, sigma_a_1, sigma_a_2);

    DXTC_TIMINGS.start("solve_all");
    if (!linear) {
      std::vector<size_t> output_dofs{2728, 3868, 4468, 929};
      solver->prepare_restricted_operator(output_dofs);
      using VectorType = typename XT::LA::Container<double, XT::LA::Backends::common_dense>::VectorType;
      auto initial_vals = solver->get_initial_values();
      RandomNumberGeneratorType rng{std::random_device()()};
      std::uniform_real_distribution<double> distribution(1, 1000);
      for (auto&& val : initial_vals)
        val *= 1e4 * distribution(rng);
      auto source_dofs = solver->restricted_op_input_dofs();
      std::cout << source_dofs << std::endl;
#if 0
    VectorType initial_vals_restr{
        3.00887845e-05, 7.40090567e-05, 7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 1.00000000e-08,
        3.79301179e-05, 1.42264780e-05, 1.51590332e-05, 1.00000000e-08, 6.04617301e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 1.00000000e-08, 3.39443780e-04, 1.51590332e-05, 1.42264780e-05, 3.79301179e-05, 1.00000000e-08,
        1.00000000e-08, 6.04617301e-05, 1.47623908e-05, 9.55564780e-06, 9.55564780e-06, 1.47623908e-05, 1.00000000e-08,
        3.30163508e-05, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.20284294e-04, 1.00000000e-08, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08,
        1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.00000000e-08, 1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02,
        3.62665787e-02, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.50692562e-04, 2.63568632e-05, 6.74968960e-05, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 2.63568632e-05, 6.74968960e-05, 1.00000000e-08,
        1.50692562e-04, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05, 7.40090567e-05, 3.39443780e-04, 1.00000000e-08};
#endif
      VectorType initial_vals_restr(source_dofs.size());
      for (size_t kk = 0; kk < source_dofs.size(); ++kk)
        initial_vals_restr[kk] = initial_vals[source_dofs[kk]];
      auto output1 = solver->apply_restricted_kinetic_operator(initial_vals_restr);
      // for (size_t ii = 0; ii < 1000; ++ii) {
      const size_t ii = 0;
      const auto actual_initial_vals = initial_vals_restr * ii / 1000.;
      output1 = solver->apply_restricted_kinetic_operator(actual_initial_vals);
      // }
      auto output2 = solver->apply_kinetic_operator(initial_vals, 0, solver->time_step_length());
      for (size_t kk = 0; kk < output_dofs.size(); ++kk)
        EXPECT_NEAR(output1[kk], output2[output_dofs[kk]], 1e-10);
    }
    const auto result = solver->solve();
    std::cout << " Result = " << std::accumulate(result.back().begin(), result.back().end(), 0.) << std::endl;
    DXTC_TIMINGS.stop("solve_all");
    parameterfile << "Elapsed time: " << DXTC_TIMINGS.walltime("solve_all") / 1000.0 << " s" << std::endl;
    parameterfile.close();

    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported: " << e.what() << std::endl;
    std::abort();
  }
} // ... main(...)
