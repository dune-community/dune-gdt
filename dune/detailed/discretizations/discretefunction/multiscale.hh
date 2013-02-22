
#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_MULTISCALE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_MULTISCALE_HH

// system
#include <sstream>
#include <vector>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

// dune-grid
#include <dune/grid/io/file/vtk/function.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/default.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/discretefunction/default.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace DiscreteFunction {

template <class MsGridImp, class LocalDiscreteFunctionImp>
class Multiscale;

template <class GridImp, class LocalDiscreteFunctionSpaceImp, class LocalVectorBackendImp>
class Multiscale<Dune::grid::Multiscale::Default<GridImp>,
                 Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst<LocalDiscreteFunctionSpaceImp,
                                                                                 LocalVectorBackendImp>>
    : public Dune::VTKFunction<typename Dune::grid::Multiscale::Default<GridImp>::GlobalGridViewType>
{
public:
  typedef Multiscale<Dune::grid::Multiscale::Default<GridImp>,
                     Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst<LocalDiscreteFunctionSpaceImp,
                                                                                     LocalVectorBackendImp>> ThisType;

  typedef Dune::grid::Multiscale::Default<GridImp> MsGridType;

  typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst<LocalDiscreteFunctionSpaceImp,
                                                                          LocalVectorBackendImp>
      LocalDiscreteFunctionType;

  //  typedef typename LocalDiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename LocalDiscreteFunctionType::ConstLocalFunctionType ConstLocalFunctionType;

  typedef typename Dune::VTKFunction<GlobalGridViewType> BaseType;

  typedef typename GlobalGridViewType::template Codim<0>::Iterator::Entity EntityType;

  typedef typename LocalDiscreteFunctionType::DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = LocalDiscreteFunctionType::dimDomain;

  static const int dimRange = LocalDiscreteFunctionType::dimRange;

  static const std::string id;

private:
  typedef typename MsGridType::EntityToSubdomainMapType EntityToSubdomainMapType;

public:
  Multiscale(const MsGridType& msGrid,
             const std::vector<Dune::shared_ptr<LocalDiscreteFunctionType>>& localDiscreteFunctions,
             const std::string name = id)
    : msGrid_(msGrid)
    , entityToSubdomainMap_(msGrid_.entityToSubdomainMap())
    , localDiscreteFunctions_(localDiscreteFunctions)
    , name_(name)
  {
    assert(localDiscreteFunctions_.size() == msGrid_.size());
  }

  const MsGridType& msGrid() const
  {
    return msGrid_;
  }

  unsigned int size() const
  {
    return msGrid_.size();
  }

  Dune::shared_ptr<const LocalDiscreteFunctionType> localDiscreteFunction(const unsigned int subdomain) const
  {
    assert(subdomain < msGrid_.size());
    return localDiscreteFunctions_[subdomain];
  }

  Dune::shared_ptr<LocalDiscreteFunctionType> localDiscreteFunction(const unsigned int subdomain)
  {
    assert(subdomain < msGrid_.size());
    return localDiscreteFunctions_[subdomain];
  }

  //  LocalFunctionType localFunction(const EntityType& entity)
  //  {
  //    // get the subdomain of this entity
  //    const unsigned int globalIndex = msGrid_.globalGridView()->indexSet().index(entity);
  //    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
  //    if (result == entityToSubdomainMap_->end()) {
  //      std::stringstream msg;
  //      msg << "Error in " << id << ": Entity " << globalIndex << " not found in the multiscale grid!";
  //      DUNE_THROW(Dune::InvalidStateException, msg.str());
  //    }
  //    const unsigned int subdomain = result->second;
  //    assert(subdomain < msGrid_.size());
  //    // get the corresponding local discrete function
  //    const LocalDiscreteFunctionType& localDiscreteFunction = *(localDiscreteFunctions_[subdomain]);
  //    // and return its local function
  //    return localDiscreteFunction.localFunction(entity);
  //  }

  ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    // get the subdomain of this entity
    const unsigned int globalIndex                                 = msGrid_.globalGridView()->indexSet().index(entity);
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": Entity " << globalIndex << " not found in the multiscale grid!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    const unsigned int subdomain = result->second;
    assert(subdomain < msGrid_.size());
    // get the corresponding local discrete function
    const LocalDiscreteFunctionType& localDiscreteFunction = *(localDiscreteFunctions_[subdomain]);
    // and return its local function
    return localDiscreteFunction.localFunction(entity);
  }

  virtual std::string name() const
  {
    return name_;
  }

  void setName(const std::string& newName = "")
  {
    name_ = newName;
  }

  /**
      @name Methods needed, to comply with the Dune::VTKFunction interface
      @{
    */
  virtual int ncomps() const
  {
    return dimRange;
  }

  virtual RangeFieldType evaluate(const int component, const EntityType& entity, const DomainType& x) const
  {
    // get the subdomain of this entity
    const unsigned int globalIndex                                 = msGrid_.globalGridView()->indexSet().index(entity);
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": Entity " << globalIndex << " not found in the multiscale grid!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    const unsigned int subdomain = result->second;
    assert(subdomain < msGrid_.size());
    // get and call the local discrete function of the entity's subdomain
    const LocalDiscreteFunctionType& localDiscreteFunction = *(localDiscreteFunctions_[subdomain]);
    return localDiscreteFunction.evaluate(component, entity, x);
  } // virtual RangeFieldType evaluate(const int component, const EntityType& entity, const DomainType& x) const
  /**
      @}
    */

private:
  const MsGridType& msGrid_;
  const Dune::shared_ptr<const EntityToSubdomainMapType> entityToSubdomainMap_;
  std::vector<Dune::shared_ptr<LocalDiscreteFunctionType>> localDiscreteFunctions_;
  std::string name_;
}; // class Multiscale

template <class GridType, class LocalDiscreteFunctionSpaceType, class LocalVectorBackendType>
const std::string
    Multiscale<Dune::grid::Multiscale::Default<GridType>,
               Dune::Detailed::Discretizations::DiscreteFunction::DefaultConst<LocalDiscreteFunctionSpaceType,
                                                                               LocalVectorBackendType>>::id =
        "detailed.discretizations.discretefunction.multiscale";

} // namespace DiscreteFunction

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_MULTISCALE_HH
