#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_EIGEN_HH

#ifdef HAVE_EIGEN

// eigen
#include <Eigen/Core>
#include <Eigen/Sparse>

// dune-common
#include <dune/common/shared_ptr.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace LA {

namespace Backend {

namespace Container {

namespace Eigen {

template <class EntryImp>
class SparseMatrix
{
public:
  typedef EntryImp EntryType;

  typedef SparseMatrix<EntryType> ThisType;

  typedef ::Eigen::SparseMatrix<EntryType> StorageType;

  SparseMatrix(unsigned int rows, unsigned int cols)
    : rows_(rows)
    , cols_(cols)
    , storage_(new StorageType(rows_, cols_))
  {
  }

  SparseMatrix(Dune::shared_ptr<StorageType> storage)
    : rows_(storage->rows())
    , cols_(storage->cols())
    , storage_(storage)
  {
  }

  SparseMatrix(const ThisType& other)
    : rows_(other.rows())
    , cols_(other.cols())
    , storage_(other.storage())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    rows_    = other.rows();
    cols_    = other.cols();
    storage_ = other.storage();
    return *this;
  }

  unsigned int rows() const
  {
    return rows_;
  }

  unsigned int cols() const
  {
    return cols_;
  }

  Dune::shared_ptr<StorageType> storage()
  {
    return storage_;
  }

  const Dune::shared_ptr<StorageType> storage() const
  {
    return storage_;
  }

  void reserve(unsigned int nnz)
  {
    storage_->reserve(nnz);
  }

  void add(unsigned int i, unsigned int j, const EntryType& val)
  {
    storage_->coeffRef(i, j) += val;
  }

  void set(unsigned int i, unsigned int j, const EntryType& val)
  {
    storage_->coeffRef(i, j) = val;
  }

  const EntryType get(unsigned int i, unsigned int j) const
  {
    return storage_->coeff(i, j);
  }

private:
  unsigned int rows_;
  unsigned int cols_;
  Dune::shared_ptr<StorageType> storage_;
}; // class SparseMatrix

template <class EntryImp>
class DenseMatrix;

template <>
class DenseMatrix<double>
{
public:
  typedef double EntryType;

  typedef DenseMatrix<EntryType> ThisType;

  typedef ::Eigen::MatrixXd StorageType;

  DenseMatrix(unsigned int rows, unsigned int cols)
    : rows_(rows)
    , cols_(cols)
    , storage_(new StorageType(rows_, cols_))
  {
  }

  DenseMatrix(Dune::shared_ptr<StorageType> storage)
    : rows_(storage->rows())
    , cols_(storage->cols())
    , storage_(storage)
  {
  }

  DenseMatrix(const ThisType& other)
    : rows_(other.rows())
    , cols_(other.cols())
    , storage_(other.storage())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    rows_    = other.rows();
    cols_    = other.cols();
    storage_ = other.storage();
    return *this;
  }

  unsigned int rows() const
  {
    return rows_;
  }

  unsigned int cols() const
  {
    return cols_;
  }

  Dune::shared_ptr<StorageType> storage()
  {
    return storage_;
  }

  const Dune::shared_ptr<StorageType> storage() const
  {
    return storage_;
  }

  void reserve()
  {
    storage_->operator=(StorageType::Zero(rows(), cols()));
  }

  void add(unsigned int i, unsigned int j, const EntryType& val)
  {
    storage_->operator()(i, j) += val;
  }

  void set(unsigned int i, unsigned int j, const EntryType& val)
  {
    storage_->operator()(i, j) = val;
  }

  const EntryType get(unsigned int i, unsigned int j) const
  {
    return storage_->operator()(i, j);
  }

private:
  unsigned int rows_;
  unsigned int cols_;
  Dune::shared_ptr<StorageType> storage_;
}; // class DenseMatrix

template <class EntryImp>
class DenseVector;

template <>
class DenseVector<double>
{
public:
  typedef double EntryType;

  typedef DenseVector<EntryType> ThisType;

  typedef ::Eigen::VectorXd StorageType;

  DenseVector(unsigned int size)
    : size_(size)
    , storage_(new StorageType(size_))
  {
  }

  DenseVector(Dune::shared_ptr<StorageType> storage)
    : size_(storage->rows())
    , storage_(storage)
  {
  }

  DenseVector(const ThisType& other)
    : size_(other.size())
    , storage_(other.storage())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    size_    = other.size();
    storage_ = other.storage();
    return *this;
  }

  unsigned int size() const
  {
    return size_;
  }

  Dune::shared_ptr<StorageType> storage()
  {
    return storage_;
  }

  const Dune::shared_ptr<StorageType> storage() const
  {
    return storage_;
  }

  void reserve()
  {
    storage_->operator=(StorageType::Zero(size()));
  }

  void add(unsigned int i, const EntryType& val)
  {
    storage_->operator()(i) += val;
  }

  void set(unsigned int i, const EntryType& val)
  {
    storage_->operator()(i) = val;
  }

  const EntryType get(unsigned int i) const
  {
    return storage_->operator()(i);
  }

private:
  unsigned int size_;
  Dune::shared_ptr<StorageType> storage_;
}; // class DenseVector

} // namespace Container

} // namespace Eigen

} // namespace Backend

} // namespace LA

} // namespace DetailedDiscretizations

} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_EIGEN_HH
