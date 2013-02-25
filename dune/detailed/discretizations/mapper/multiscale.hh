
#ifndef DUNE_DETAILEL_DISCRETIZATIONS_MAPPER_MULTISCALE_HH
#define DUNE_DETAILEL_DISCRETIZATIONS_MAPPER_MULTISCALE_HH

// system
#include <map>
#include <sstream>
#include <vector>

// dune-common
#include <dune/common/exceptions.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Mapper {

template <class IndexImp = size_t>
class Multiscale
{
public:
  typedef IndexImp IndexType;

  typedef Multiscale<IndexType> ThisType;

  static const std::string id;

  Multiscale()
    : subdomainMap_(0)
    , numSubdomains_(0)
    , size_(0)
    , prepared_(false)
    , finalized_(false)
  {
  }

  void prepare()
  {
    assert(!finalized_);
    if (!prepared_) {
      subdomainMap_ = new std::map<unsigned int, IndexType>();
      prepared_     = true;
    }
    return;
  } // void prepare()

  void add(const unsigned int subdomain, const IndexType size)
  {
    assert(prepared_);
    assert(!finalized_);
    assert(0 <= subdomain);
    assert(0 <= size);
    // create entry for this subdomain if needed (doing this explicitly instead of just using insert only to increment
    // size)
    if (subdomainMap_->find(subdomain) == subdomainMap_->end()) {
      subdomainMap_->insert(std::pair<unsigned int, IndexType>(subdomain, size));
      ++numSubdomains_;
    }
    return;
  } // void add(const unsigned int subdomain, const IndexType size)

  void finalize()
  {
    assert(prepared_);
    if (!finalized_) {
      // test for consecutive numbering of subdomains and create vector of local sizes
      for (unsigned int subdomain = 0; subdomain < numSubdomains_; ++subdomain) {
        const typename std::map<unsigned int, IndexType>::const_iterator result = subdomainMap_->find(subdomain);
        if (result == subdomainMap_->end()) {
          std::stringstream msg;
          msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
          DUNE_THROW(Dune::InvalidStateException, msg.str());
        } else {
          const IndexType& localSize = result->second;
          localSizes_.push_back(localSize);
          globalStartIndices_.push_back(size_);
          size_ += localSize;
        }
      } // // test for consecutive numbering of subdomains and create vector of local sizes
      // clean up
      delete subdomainMap_;
      finalized_ = true;
    }
    return;
  } // void finalize()

  IndexType numSubdomains() const
  {
    assert(finalized_);
    return numSubdomains_;
  }

  IndexType size() const
  {
    assert(finalized_);
    return size_;
  }

  IndexType toGlobal(const unsigned int subdomain, const IndexType localIndex) const
  {
    assert(finalized_);
    assert(subdomain < numSubdomains());
    assert(0 <= localIndex);
    assert(localIndex < localSizes_[subdomain]);
    return globalStartIndices_[subdomain] + localIndex;
  } // IndexType toGlobal(const unsigned int subdomain, const InedxType localIndex) const

private:
  std::map<unsigned int, IndexType>* subdomainMap_;
  unsigned int numSubdomains_;
  unsigned int size_;
  std::vector<IndexType> localSizes_;
  std::vector<IndexType> globalStartIndices_;
  bool prepared_;
  bool finalized_;
}; // class Multiscale

template <class IndexType>
const std::string Multiscale<IndexType>::id = "detailed.discretizations.mapper.multiscale";

} // namespace Mapper

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETAILEL_DISCRETIZATIONS_MAPPER_MULTISCALE_HH
