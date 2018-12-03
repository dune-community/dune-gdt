// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DOF_VECTOR_HH
#define DUNE_GDT_DISCRETEFUNCTION_DOF_VECTOR_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/local/dof-vector.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>

namespace Dune {
namespace GDT {


template <class Vector, class GridView>
class ConstDofVector
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

  using ThisType = ConstDofVector<Vector, GridView>;

public:
  using VectorType = Vector;
  using ConstLocalDofVectorType = ConstLocalDofVector<Vector, GridView>;
  using MapperType = MapperInterface<GridView>;
  using ElementType = typename MapperType::ElementType;

  ConstDofVector(const MapperType& mapper, const VectorType& vec)
    : mapper_(mapper)
    , vector_(vec)
  {
    if (vector_.size() != mapper_.size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "mapper_.size() = " << mapper_.size() << "\n   "
                                     << "vector_.size() = " << vector_.size());
  }

  const VectorType& vector() const
  {
    return vector_;
  }

  ConstLocalDofVectorType localize() const
  {
    return ConstLocalDofVectorType(mapper_, vector_);
  }

protected:
  const MapperType& mapper_;

private:
  const VectorType& vector_;
}; // class ConstDofVector


template <class Vector, class GridView>
class DofVector : public ConstDofVector<Vector, GridView>
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

  using ThisType = DofVector<Vector, GridView>;
  using BaseType = ConstDofVector<Vector, GridView>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::MapperType;
  using typename BaseType::VectorType;
  using LocalDofVectorType = LocalDofVector<Vector, GridView>;

  DofVector(const MapperType& mapper, VectorType& vec)
    : BaseType(mapper, vec)
    , vector_(vec)
  {}

  using BaseType::vector;

  VectorType& vector()
  {
    return vector_;
  }

  using BaseType::localize;

  LocalDofVectorType localize()
  {
    return LocalDofVectorType(mapper_, vector_);
  }

  /// \name These methods are required for grid adaptation.
  /// \{

  void resize_after_adapt()
  {
    vector_.resize(mapper_.size());
  }

  /// \}

private:
  using BaseType::mapper_;
  VectorType& vector_;
}; // class DofVector


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DOF_VECTOR_HH
