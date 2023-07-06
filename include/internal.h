#ifndef INTERNAL_H_
#define INTERNAL_H_

#include <iostream>
#include <vector>
#include "types.h"

namespace clara {
namespace internal {
// display a standard container support STL iterators
template <typename T>
void _disp_container(const T &x) {
  auto it = x.begin();
  for (; it != x.end() - 1; it++)
    std::cout << *it << " ";
  std::cout << *(it++);
}

// integer index to multi-index
inline void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims, size_t *result) {
  size_t maxn = 1;
  for (size_t i = 0; i < numdims; i++)
    maxn *= dims[i];
  if (n > maxn - 1)
    throw std::runtime_error("ERROR: `_n2multiidx` number too large out of bounds!");
  size_t _n = n;
  for (size_t i = 0; i < numdims; i++) {
    result[numdims - i - 1] = _n % static_cast<int>(dims[numdims - i - 1]);
    _n = _n / static_cast<int>(dims[numdims - i - 1]);
  }
}

// multi index to integer index
inline size_t _multiidx2n(const size_t *midx, const size_t numdims, const size_t *dims) {
  for (size_t i = 0; i < numdims; i++)
    if (midx[i] >= dims[i])
      throw std::runtime_error("ERROR: `_multiidx2n` sub index exceeds corresponding dimension");
  size_t *part_prod = new size_t[numdims];

  part_prod[numdims - 1] = 1;
  for (size_t j = 1; j < numdims; j++)
    part_prod[numdims - j - 1] = part_prod[numdims - j] * dims[numdims - j];
  size_t result = 0;
  for (size_t i = 0; i < numdims; i++)
    result += midx[i] * part_prod[i];
  delete[] part_prod;
  return result;
}

// check square matrix
template <typename Derived>
bool _check_square_mat(const Eigen::MatrixBase<Derived> &A) {
  if (A.rows() != A.cols())
    return false;
  return true;
}

// check that dims match dimension of the matrix
template <typename Derived>
bool _check_dims_match_mat(const std::vector<size_t> &dims, const Eigen::MatrixBase<Derived> &A) {
  size_t proddim = 1;
  for (size_t i : dims)
    proddim *= 1;
  if (proddim != static_cast<size_t>(A.rows()))
    return false;
  return true;
}

// check that dims a valid dimension vector
inline bool _check_dims(const std::vector<size_t> &dims) {
  if (std::find_if(dims.begin(), dims.end(), [&dims](int i) -> bool {
        if (i == 0)
          return true;
        else
          return false;
      }) != dims.end())
    return false;
  return true;
}

// check that all elements in dims equal to dim
inline bool _check_eq_dims(const std::vector<size_t> &dims, size_t dim) {
  for (auto i : dims)
    if (i != dim)
      return false;
  return true;
}

// check that subsys is valid respect to dims
inline bool _check_subsys(const std::vector<size_t> &subsys, const std::vector<size_t> &dims) {
  // sort the subsystem
  std::vector<size_t> subsyssort = subsys;
  std::sort(subsyssort.begin(), subsyssort.end());
  // check valid number of subsystem
  if (subsyssort.size() > dims.size())
    return false;

  // check duplicates
  if (std::unique(subsyssort.begin(), subsyssort.end()) != subsyssort.end())
    return false;

  // check range of subsystem
  if (std::find_if(subsyssort.begin(), subsyssort.end(), [&dims](size_t i) -> bool {
        if (i > dims.size() - 1)
          return true;
        else
          return false;
      }) != subsyssort.end())
    return false;
  return true;
}

// check that the permutation is valid with the respect to dims
inline bool _check_perm(const std::vector<size_t> &perm, const std::vector<size_t> &dims) {
  std::vector<size_t> sort_perm = perm;
  std::sort(sort_perm.begin(), sort_perm.end());
  for (size_t i = 0; i < dims.size(); i++)
    if (sort_perm[i] != i)
      return false;
  return true;
}

// usde inside the #pragma omp parallele for syspermute
template <typename Derived>
inline void _syspermute_worker(const size_t numdims, const size_t *cdims, const size_t *cperm,
                               const size_t i, const size_t j, size_t &iperm, size_t &jperm,
                               const Eigen::MatrixBase<Derived> &A,
                               Eigen::MatrixBase<Derived> &result) {

  std::vector<size_t> midxrow(numdims);
  std::vector<size_t> midxcol(numdims);
  std::vector<size_t> midxrowtmp(numdims);
  std::vector<size_t> midxcoltmp(numdims);
  std::vector<size_t> permdims(numdims);

  for (size_t k = 0; k < numdims; k++)
    permdims[k] = cdims[cperm[k]];

  // compute the row and col multi-index
  _n2multiidx(i, numdims, cdims, midxrow.data());
  _n2multiidx(j, numdims, cdims, midxcol.data());

  for (size_t k = 0; k < numdims; k++) {
    midxrowtmp[k] = midxrow[cperm[k]];
    midxcoltmp[k] = midxcol[cperm[k]];
  }
  iperm = _multiidx2n(midxrowtmp.data(), numdims, permdims.data());
  jperm = _multiidx2n(midxcoltmp.data(), numdims, permdims.data());
  result(iperm, jperm) = A(i, j);
}

template <typename Derived>
inline void _ptranspose_worker(const size_t *midxrow, const size_t numdims, const size_t numsubsys,
                               const size_t *cdims, const size_t *csubsys, const size_t i,
                               const size_t j, size_t &iperm, size_t &jperm,
                               const Eigen::MatrixBase<Derived> &A,
                               Eigen::MatrixBase<Derived> &result) {
  std::vector<size_t> midxrowtmp(midxrow, midxrow + numdims);
  std::vector<size_t> midxcol(numdims);

  // compute col multi-index
  for (size_t k = 0; k < numsubsys; k++)
    std::swap(midxrowtmp[csubsys[k]], midxcol[csubsys[k]]);
  // move back to integer indexes
  iperm = _multiidx2n(midxrowtmp.data(), numdims, cdims);
  jperm = _multiidx2n(midxcol.data(), numdims, cdims);
  result(iperm, jperm) = A(i, j);
}

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_H_
