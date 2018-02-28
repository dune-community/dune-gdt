template <class GV>
class BasisVisualization : public Dune::VTKFunction<GV>
{
  static_assert(dimRangeCols == 1, "Not implemented for matrixvalued spaces yet!");

public:
  typedef typename GlobalBasisType::DomainType DomainType;
  typedef typename GlobalBasisType::RangeType RangeType;

  BasisVisualization(const derived_type& sp, const size_t ind, const std::string nm = "basis")
    : space_(sp)
    , index_(ind)
    , values_(space_.mapper().maxNumDofs(), RangeType(0))
    , name_(nm)
  {
    if (index_ >= space_.mapper().maxNumDofs())
      DUNE_THROW(Dune::RangeError,
                 "index has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << index_ << ")!");
  }

  virtual std::string name() const
  {
    return name_;
  }

  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  /* @{ */
  virtual int ncomps() const
  {
    return dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& xx) const
  {
    const auto baseFunctionSet = space_.base_function_set(entity);
    if (component < 0)
      DUNE_THROW(Dune::RangeError, "component must not be negative (is " << component << ")!");
    if (component < boost::numeric_cast<int>(baseFunctionSet.size())) {
      baseFunctionSet.evaluate(xx, values_);
      assert(component < boost::numeric_cast<int>(values_.size()) && "This should not happen!");
      return values_[index_][component];
    } else if (component < boost::numeric_cast<int>(space_.mapper().maxNumDofs()))
      return 0.0;
    else
      DUNE_THROW(Dune::RangeError,
                 "component has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << component << ")!");
  }
  /* @} */

private:
  const derived_type& space_;
  const size_t index_;
  mutable std::vector<RangeType> values_;
  const std::string name_;
}; // class BasisVisualization


void visualize(const std::string filename_prefix = "") const
{
  const std::string filename = filename_prefix.empty() ? "dune.gdt.space" : filename_prefix;
  const auto tmp_storage = XT::Grid::make_tmp_view(grid_layer());
  const auto& grd_vw = tmp_storage.access();
  using GridViewType = std::decay_t<decltype(grd_vw)>;
  VTKWriter<GridViewType> vtk_writer(grd_vw, Dune::VTK::nonconforming);
  for (size_t ii = 0; ii < mapper().maxNumDofs(); ++ii) {
    std::string number = "";
    if (ii == 1)
      number = "1st";
    else if (ii == 2)
      number = "2nd";
    else if (ii == 3)
      number = "3rd";
    else
      number = XT::Common::to_string(ii) + "th";
    const auto iith_baseFunction =
        std::make_shared<BasisVisualization<GridViewType>>(this->as_imp(*this), ii, number + " basis");
    vtk_writer.addVertexData(iith_baseFunction);
  }
  vtk_writer.write(filename);
} // ... visualize(...)
