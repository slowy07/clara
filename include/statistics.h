#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <type_traits>
#include <vector>

#include "classFunction/exception.h"
#include "internal/util.h"
#include "traits.h"
#include "types.h"
namespace clara {

/*
 * @brief uniform probability distribution vector
 * @param N size of the alphabet
 * @return Real vector consiting of a uniform distribution of size N
 */
inline std::vector<double> uniform(idx N) {
  if (N == 0)
    throw Exception("clara::uniform", Exception::Type::ZERO_SIZE);
  return std::vector<double>(N, 1. / N);
}

/**
 * @brief marginal distribution
 * @param probXY real matrix representing the joint probability distribution
 * of x nd y in lexicographical order
 */
inline std::vector<double> marginalX(const dmat& probXY) {
  if (!internal::check_nonzero_size(probXY))
    throw Exception("clara::marginalX", Exception::Type::ZERO_SIZE);

  std::vector<double> result(probXY.rows(), 0);
  for (idx i = 0; i < static_cast<idx>(probXY.rows()); ++i) {
    for (idx j = 0; j < static_cast<idx>(probXY.cols()); ++j) {
      result[i] += probXY(i, j);
    }
  }
  return result;
}

/**
 * @brief marginal distribution
 * @param probXY real matrix representing the joint probability distribution
 * of X and Y in lexicographical order
 * @return real vector consiting of the marginal distribution of Y
 */
inline std::vector<double> marginalY(const dmat& probXY) {
  if (!internal::check_nonzero_size(probXY))
    throw Exception("clara::marginalY", Exception::Type::ZERO_SIZE);
  return marginalX(probXY.transpose());
}

/**
 * @brief average
 * @param prob real probability vector representing the probability distribtuib of X
 * @param X random variable values repreented by an STL-like container
 * @return average of X
 */
template <typename Container>
double avg(const std::vector<double>& prob, const Container& X,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::avg", Exception::Type::ZERO_SIZE);
  if (!internal::check_matching_sizes(prob, X))
    throw Exception("clara::avg", Exception::Type::SIZE_MISMATCH);

  double result = 0;
  for (idx i = 0; i < prob.size(); ++i)
    result += prob[i] * X[i];
  return result;
}

/**
 * @brief covariance
 * @param probXY real matrix representing the join probability distribution
 * of X and Y in lexicographical order (X labels the rows Y labels the column)
 * @param X random variable values represented by an STL-like container
 * @param Y random variable values represented by an STL-like container
 * @return covariance of X and Y
 */
template <typename Container>
double cov(const dmat& probXY, const Container& X, const Container& Y,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(X) || !internal::check_nonzero_size(Y))
    throw Exception("clara::cov", Exception::Type::ZERO_SIZE);
  if (static_cast<idx>(probXY.rows()) != X.size() || static_cast<idx>(probXY.cols()) != Y.size())
    throw Exception("clara::cov", Exception::Type::SIZE_MISMATCH);
  std::vector<double> probX = marginalX(probXY);
  std::vector<double> probY = marginalY(probXY);

  double result = 0;
  for (idx i = 0; i < X.size(); ++i) {
    for (idx j = 0; j < Y.size(); ++j) {
      result += probXY(i, j) * X[i] * Y[j];
    }
  }
  return result - avg(probX, X) * avg(probY, Y);
}

template <typename Container>
double var(const std::vector<double>& prob, const Container& X,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::var", Exception::Type::ZERO_SIZE);
  if (!internal::check_matching_sizes(prob, X))
    throw Exception("clara::var", Exception::Type::SIZE_MISMATCH);
  Eigen::VectorXcd diag(prob.size());
  for (idx i = 0; i < prob.size(); ++i)
    diag(i) = prob[i];
  dmat probXX = diag.asDiagonal();
  return cov(probXX, X, X);
}

template <typename Container>
double sigma(const std::vector<double>& prob, const Container& X,
             typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw Exception("clara::sigma", Exception::Type::ZERO_SIZE);
  if (!internal::check_matching_sizes(prob, X))
    throw Exception("clara::sigma", Exception::Type::SIZE_MISMATCH);
  return std::sqrt(var(prob, X));
}

template <typename Container>
double cor(const dmat& probXY, const Container& X, const Container& Y,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(X) || !internal::check_nonzero_size(Y))
    throw Exception("clara::cor", Exception::Type::ZERO_SIZE);
  if (static_cast<idx>(probXY.rows()) != X.size() || static_cast<idx>(probXY.cols()) != Y.size())
    throw Exception("clara::cor", Exception::Type::SIZE_MISMATCH);
  return cov(probXY, X, Y) / (sigma(marginalX(probXY), X) * sigma(marginalX(probXY), Y));
}

}  // namespace clara

#endif  // !STATISTICS_H_
