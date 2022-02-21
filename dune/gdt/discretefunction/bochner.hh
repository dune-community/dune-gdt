// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_BOCHNER_HH
#define DUNE_GDT_DISCRETEFUNCTION_BOCHNER_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/la/container/vector-array/list.hh>
#include <dune/xt/grid/search.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/bochner.hh>

namespace Dune {
namespace GDT {


/**
 * \todo Turn this into a parametric LocalizableFunction.
 * \todo add ConstDiscreteBochnerFunction
 */
template <class V, class GV, size_t r = 1, size_t rC = 1, class R = double>
class DiscreteBochnerFunction
{
public:
  DiscreteBochnerFunction(const BochnerSpace<GV, r, rC, R>& bochner_space,
                          XT::LA::ListVectorArray<V>& dof_vectors,
                          const std::string nm = "")
    : bochner_space_(bochner_space)
    , dof_vectors_(dof_vectors)
    , name_(nm.empty() ? "DiscreteBochnerFunction" : nm)
  {
    DUNE_THROW_IF(this->dof_vectors().length() != bochner_space_.temporal_space().mapper().size(),
                  Exceptions::space_error,
                  "\n   this->dof_vectors().length() = " << this->dof_vectors().length() << "\n   "
                                                         << bochner_space_.temporal_space().mapper().size());
    for (const auto& vec : this->dof_vectors())
      DUNE_THROW_IF(vec.vector().size() != bochner_space_.spatial_space().mapper().size(),
                    Exceptions::space_error,
                    "\n   vec.vector().size() = " << vec.vector().size() << "\n   "
                                                  << bochner_space_.spatial_space().mapper().size());
  } // DiscreteBochnerFunction(...)

  DiscreteBochnerFunction(const BochnerSpace<GV, r, rC, R>& bochner_space,
                          XT::LA::ListVectorArray<V>*&& dof_vectors,
                          const std::string nm = "")
    : bochner_space_(bochner_space)
    , dof_vectors_(std::move(dof_vectors))
    , name_(nm.empty() ? "DiscreteBochnerFunction" : nm)
  {
    DUNE_THROW_IF(this->dof_vectors().length() != bochner_space_.temporal_space().mapper().size(),
                  Exceptions::space_error,
                  "\n   this->dof_vectors().length() = " << this->dof_vectors().length() << "\n   "
                                                         << bochner_space_.temporal_space().mapper().size());
    for (const auto& vec : this->dof_vectors())
      DUNE_THROW_IF(vec.vector().size() != bochner_space_.spatial_space().mapper().size(),
                    Exceptions::space_error,
                    "\n   vec.vector().size() = " << vec.vector().size() << "\n   "
                                                  << bochner_space_.spatial_space().mapper().size());
  } // DiscreteBochnerFunction(...)

  DiscreteBochnerFunction(const BochnerSpace<GV, r, rC, R>& bochner_space, const std::string nm = "")
    : bochner_space_(bochner_space)
    , dof_vectors_(new XT::LA::ListVectorArray<V>(bochner_space_.spatial_space().mapper().size(),
                                                  bochner_space_.temporal_space().mapper().size()))
    , name_(nm.empty() ? "DiscreteBochnerFunction" : nm)
  {}

  const BochnerSpace<GV, r, rC, R>& space() const
  {
    return bochner_space_;
  }

  const XT::LA::ListVectorArray<V>& dof_vectors() const
  {
    return dof_vectors_.access();
  }

  XT::LA::ListVectorArray<V>& dof_vectors()
  {
    return dof_vectors_.access();
  }

  std::string name() const
  {
    return name_;
  }

  DiscreteFunction<V, GV, r, rC, R> evaluate(const double& time) const
  {
    const double t =
        XT::Common::clamp(time, bochner_space_.time_interval().first, bochner_space_.time_interval().second);
    const auto search_result =
        XT::Grid::make_entity_in_level_search(bochner_space_.temporal_space().grid_view())(std::vector<double>(1, t));
    DUNE_THROW_IF(!search_result.at(0),
                  Exceptions::finite_element_error,
                  "when evaluating timedependent function: could not find time_interval for time = " << t);
    const auto& time_interval = *search_result.at(0);
    const auto temporal_basis = bochner_space_.temporal_space().basis().localize(time_interval);
    const auto time_in_reference_element = time_interval.geometry().local(t);
    const auto temporal_basis_values = temporal_basis->evaluate_set(time_in_reference_element);
    V result(bochner_space_.spatial_space().mapper().size(), 0.);
    const auto global_dof_indices = bochner_space_.temporal_space().mapper().global_indices(time_interval);
    for (size_t ii = 0; ii < temporal_basis->size(); ++ii)
      result.axpy(temporal_basis_values[ii], this->dof_vectors()[global_dof_indices[ii]].vector());
    return make_discrete_function(bochner_space_.spatial_space(), std::move(result));
  } // ... evaluate(...)

  void visualize(const std::string filename_prefix, const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    DUNE_THROW_IF(
        filename_prefix.empty(), XT::Common::Exceptions::wrong_input_given, "filename_prefix must not be empty!");
    bool use_counter = false;
    size_t counter = 0;
    for (const auto& annotated_vector : dof_vectors_.access())
      if (!annotated_vector.note().has_key("_t"))
        use_counter = true;
    for (const auto& annotated_vector : dof_vectors_.access()) {
      auto df = make_discrete_function(bochner_space_.spatial_space(), annotated_vector.vector(), name_);
      if (use_counter) {
        df.visualize(filename_prefix + "_" + XT::Common::to_string(counter), vtk_output_type);
        ++counter;
      } else {
        const double time = annotated_vector.note().get("_t").at(0);
        df.visualize(filename_prefix + "_" + XT::Common::to_string(time), vtk_output_type);
      }
    }
  } // ... visualize(...)

private:
  const BochnerSpace<GV, r, rC, R>& bochner_space_;
  XT::Common::StorageProvider<XT::LA::ListVectorArray<V>> dof_vectors_;
  const std::string name_;
}; // class DiscreteBochnerFunction


template <class GV, size_t r, size_t rC, class R, class V>
DiscreteBochnerFunction<V, GV, r, rC, R> make_discrete_bochner_function(const BochnerSpace<GV, r, rC, R>& bochner_space,
                                                                        XT::LA::ListVectorArray<V>& dof_vectors,
                                                                        const std::string name = "")
{
  return DiscreteBochnerFunction<V, GV, r, rC, R>(bochner_space, dof_vectors, name);
}


template <class VectorType, class GV, size_t r, size_t rC, class R>
typename std::enable_if<XT::LA::is_vector<VectorType>::value, DiscreteBochnerFunction<VectorType, GV, r, rC, R>>::type
make_discrete_bochner_function(const BochnerSpace<GV, r, rC, R>& bochner_space, const std::string name = "")
{
  return DiscreteBochnerFunction<VectorType, GV, r, rC, R>(bochner_space, name);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_BOCHNER_HH
