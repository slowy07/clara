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
 * @brief generate a uniform probability distribution vector
 *
 * this function generate a real vector of size N representing a uniform probability distribution.
 * in a uniform probability distribution, all elements have an equal probability of occurence.
 * each elements in the generated vector will have a value of 1/N, where N is the size of the
 * alphabet the function returns  a vector where each element is set to 1/N, ensuring that the
 * vector representing a uniform distribution
 *
 * @param N the size of the alphabet for which the uniform distribution is generated
 * @return real vector representing a uniform probability distribution of size N
 *
 * @throws exception::ZeroSize if 'N' is zero, indicating an invalid vector size
 *
 * @example
 * // usage of the uniform function to generate a uniform probability distribution vector of size
 * int N = 4;
 * std::vector<double> uniform_dist = uniform(N);
 */
inline std::vector<double> uniform(idx N) {
  if (N == 0)
    throw exception::ZeroSize("clara::uniform()");
  return std::vector<double>(N, 1. / N);
}

/**
 * @brief compute the marginal distribution for a given joint probability distribution
 *
 * this function calculates the marginal distribution of the variable 'X' forma given joint
 * probability distribution. the input 'probXY' is real matrix representing the joint probability
 * distribution of two matrices 'X' and 'Y' in lexicographical order, where 'probXy(i, j)' represent
 * the probability of 'X=i' and 'Y=j'
 *
 * @param probXY a real matrix representing the joint probability distribution of 'X' and 'Y'
 * @return real vector representing the matrignal disribution of 'X'
 *
 * @throws exception::ZeroSize if 'probXY' is an empty matrix, indicating invalid input
 *
 * @example
 * // usage of marginalX function to calculate the marginal distribution of variable 'X'
 * dmat probXY = ...;
 * std::vector<double> marginalDistX = marginalX(probXY);
 */
inline std::vector<double> marginalX(const dmat& probXY) {
  if (!internal::check_nonzero_size(probXY))
    throw exception::ZeroSize("clara::marginalX()");

  std::vector<double> result(probXY.rows(), 0);
  for (idx i = 0; i < static_cast<idx>(probXY.rows()); ++i) {
    for (idx j = 0; j < static_cast<idx>(probXY.cols()); ++j) {
      result[i] += probXY(i, j);
    }
  }
  return result;
}

/**
 * @brief compute the marginal distribution of variable 'Y' from a given joint probability
 * distribution
 *
 * this function calculates the marginal distribution of the variable 'Y' from a given joint
 * probability distribution. the input 'probXY' is a real matrix representing the joint probability
 * distribution of two variables 'X' and 'Y' in lexicographical order, where 'probXY(i, j)'
 * represents the probability of 'X=i' and 'Y=j' the function function returns a real vector
 * containing the marginal distribution of 'Y', where each element 'j' represents the probability of
 * 'Y=j'. the marginal distribution of 'Y' is obtained by summing up the probabilities of all events
 * where 'Y=j' across all possible values of 'X'
 *
 * @param probXY real matrix representing the joint probability distribution of 'X' and 'Y'
 * @return real vector representing the marginal distribution of 'Y'
 * @throws exception::ZeroSize 'probXY' is an empty matrix, indicating an invalid input
 *
 * @example
 * // usage of marginalY functio to calculate the marginal distribution of variable 'Y'
 * dmat probXY = ...;
 * std::vector<double> marginalDistY = marginalY(probXY);
 */
inline std::vector<double> marginalY(const dmat& probXY) {
  if (!internal::check_nonzero_size(probXY))
    throw exception::ZeroSize("clara::marginalY()");
  return marginalX(probXY.transpose());
}

/**
 * @brief calculate the averange of a random variable 'X' based on its probability distribution
 * this function calculate the averange (or expected value) of a random variable 'X' based on its
 * probability distribution the averange 'X' is given the sum of the products of the values of 'X'
 * and their corresponding probabilities the probabilities are provided as a real probability vector
 * 'prob', and the values of 'X' are represented by an STL-like container 'X'
 *
 * @param prob real probability vector representing the probability distribution of 'X'
 * @param X an STL-like container representing the value of the random variable 'X'
 *
 * @return the average of the random variable 'X' based on its probability distribution
 *
 * @throws exception::ZeroSize if 'prob' or 'X' is an empty container, indicating an ivalid input
 * @throws exception::SizeMismatch if 'prob' and 'X' have different sizes, indicating an invalid
 * input
 *
 * @example
 * // usage of avg function to calculate the average of a random variable 'X'
 * std::vector<double> prob = {0.2, 0.3, 0.5};
 * std::vector<double> x = {1.0, 2.0, 3.0};
 * double averageX = avg(prob, X);
 */
template <typename Container>
double avg(const std::vector<double>& prob, const Container& X,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::avg()");
  if (!internal::check_matching_sizes(prob, X))
    throw exception::SizeMismatch("clara::avg()");

  double result = 0;
  for (idx i = 0; i < prob.size(); ++i)
    result += prob[i] * X[i];
  return result;
}

/**
 * @brief calculate the covariance between two random variables 'X' and 'Y' based on their join
 * probability distribution
 *
 * this function calculate the covariance between two random variables 'X' and 'Y' using their joint
 * probability distribution 'probXY' joint probability distribution is represented as a real matrix
 * 'probXY', where 'probXY(i, j)' represents the joint probability of 'X[i]' and 'Y[j]'
 *
 * @tparam Container an STL-like container type representing the values of 'X' and 'Y'
 *
 * @param probXY real matrix representing the joint probability distribution of 'X' and 'Y' in
 * lexicographical order  'probXY(i, j)' represents the joint probability of 'X[i]' and 'Y[j]'
 * @param X random variable values represented by an STL-like container
 * @param Y random variable values represented by an STL-like container
 *
 * @return the covariance between 'X' and 'Y' based on their joint probability distribution
 *
 * @throws exception::ZeroSize if either 'X' or 'Y' has zero size, indicating an invalid input size
 * @throws exception::SizeMismatch if the size of 'probXY' does not match the size of 'X' and 'Y',
 *                                 indicating an inconsistency in the data
 *
 * @example
 * // usage of the cov function to calculate the covariance between two random variables 'X' and 'Y'
 * idx D = 3;
 * std::vector<double> probVector = {0.1, 0.2, 0.3};
 * dmat probMatrix = dmat::Zero(D, D);
 * probMatrix(0, 1) = 0.2;
 * probMatrix(1, 2) = 0.3;
 * probMatrix(2, 0) = 0.1;
 *
 * std::vector<double> X = {1, 2, 3};
 * std::vector<double> Y = {4, 5, 6};
 * double covarianceXY = cov(probMatrix, X, Y);
 */
template <typename Container>
double cov(const dmat& probXY, const Container& X, const Container& Y,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(X) || !internal::check_nonzero_size(Y))
    throw exception::ZeroSize("clara::cov()");
  if (static_cast<idx>(probXY.rows()) != X.size() || static_cast<idx>(probXY.cols()) != Y.size())
    throw exception::SizeMismatch("clara::cov()");
  // calculate thr marginal probability distribution of 'X' and 'Y'
  std::vector<double> probX = marginalX(probXY);
  std::vector<double> probY = marginalY(probXY);

  // calculate the covarience between 'X' and 'Y'
  double result = 0;
  for (idx i = 0; i < X.size(); ++i) {
    for (idx j = 0; j < Y.size(); ++j) {
      result += probXY(i, j) * X[i] * Y[j];
    }
  }
  return result - avg(probX, X) * avg(probY, Y);
}

/**
 * @brief calculate the variance of random variable 'X' based on its probability distribution
 *
 * this function calculates the veariance of a random variable 'X' using its probability
 * distribution 'prob' the probability distribution of 'X' is represented as a real vector 'prob',
 * where 'prob[i]' represent the probability of event 'X = x[i]'
 *
 * @tparam Container an STL-like container type representing the value of 'X'
 * @param prob A real vector representing the probability distribution of 'X'. 'prob[i]' represente
 * the probability of event 'X = X[i]'
 * @param X random variable values represented by an STL-like container
 * @return the variance of 'X' based on its probability distribution
 *
 * @throws exception::ZeroSize if either prob or x has zero size, indicating an invalid input size
 * @throws exception::SizeMismatch if the size of 'prob' does not match the size of 'X' indicating
 * an inconsistency in the data
 *
 * @example
 * // usage of the var function to calculate the variance of a random variable 'X'
 * std::vector<double> prob = {0.1, 0.2, 0.3, 0.4};
 * std::vector<double> X = {1, 2, 3, 4};
 * double varianceX = var(prob, x);
 */
template <typename Container>
double var(const std::vector<double>& prob, const Container& X,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::var()");
  if (!internal::check_matching_sizes(prob, X))
    throw exception::SizeMismatch("clara::var()");

  // create an diagonal matrix with 'prob' as the diagonal elements
  Eigen::VectorXcd diag(prob.size());
  for (idx i = 0; i < prob.size(); ++i)
    diag(i) = prob[i];
  dmat probXX = diag.asDiagonal();

  // calculate the variance of 'X' using the covariance function 'cov'
  return cov(probXX, X, X);
}

/**
 * @brief calculate the standard deviation (sigma) of a random variable 'X' based on its probability
 * distribution
 *
 * this function calculates the standard deviation of a random variable 'X' using its probability
 * distribution 'rpob' the probability distribution of 'X' is represented as a real vector 'prob',
 * where 'prob[i]' represente the probability of the event 'X = X[i]', the standard deviation
 * (sigma).
 *
 * @tparam Container an STL-like container type representing the value of 'X'
 * @param prob real vector representing the probability distribution of 'X'. 'prob[i]' repersente
 * the probability of event 'X = X[i]'
 * @param X random variable values represented by an STL-like container
 *
 * @return the standard deviation (sigma) of 'X' basedd on its probability distribution
 *
 * @throws exception::ZeroSize if either 'prob' or 'X' has zero size, indicating an invalid input
 * size.
 * @throws exception::SizeMismatch if the size of 'prob' does not match the size of 'X', indicating
 * an inconsistency in the data
 *
 * @example
 * // usage the sigma function to calculate the standard deviation of random variable 'X'
 * std::vector<double> prob {0.1, 0.2, 0.3, 0.4};
 * std::vector<double> X = {1, 2, 3, 4}
 * double standardDeviationX = sigma(prob, X);
 */
template <typename Container>
double sigma(const std::vector<double>& prob, const Container& X,
             typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(prob))
    throw exception::ZeroSize("clara::sigma()");
  if (!internal::check_matching_sizes(prob, X))
    throw exception::SizeMismatch("clara::sigma()");
  // calculate the standard deviation of 'X' using the square root of its variance
  return std::sqrt(var(prob, X));
}

/**
 * @brief calculate the correlation coefficient between two random variable 'X' amd 'Y' based on
 * their joint probability distribution this function calculates the correlation coefficient between
 * two random variables 'X' and 'Y' using their joint probability distribution 'probXY', joint
 * probability distribution of 'X' and 'Y' is represented as a real matrix 'probXY', where
 * 'probXY(i, j)' represents the probability of the event 'X = X[i]' and 'Y = Y[j]'
 *
 * @tparam Container an STL-like containertype representing the values of 'X' and 'Y'
 * @param probXY real matrix representing the joint probability distribution of 'X' and 'Y' in
 * lexicographical order, 'probXY(i, j)' represents the probability of the event 'X = X[i]' and 'Y =
 * Y[j]'
 * @param X random variable values represented by an STL-like container
 * @param Y random variable values represented by an STL-like container
 *
 * @throws exception::ZeroSize if either 'X' or 'Y' has zero size, indicating an invalid input size
 * @throws exception::SizeMismatch  if the size of 'probXY' does not match the sizes of 'X' and 'Y',
 * indicating an inconsistency the data
 *
 * @example
 * // usage the cor function to calculate the correlation coefficient between two random variables
 * 'X' and 'Y' std::vector<double> prob = {0.1, 0.2, 0.3, 0.4};
 * std::vector<double> X = {1, 2, 3, 4};
 * std::vector<double> Y = {2, 4, 6, 8};
 * dmat probXY = { {0.05, 0.10, 0.15, 0.20}, {0.10, 0.20,
 * 0.30, 0.40}, {0.15, 0.30, 0.45, 0.60}, {0.20, 0.40, 0.60, 0.80} };
 *
 * double correlationXY = cor(probXY, X, Y);
 */
template <typename Container>
double cor(const dmat& probXY, const Container& X, const Container& Y,
           typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  if (!internal::check_nonzero_size(X) || !internal::check_nonzero_size(Y))
    throw exception::ZeroSize("clara::cor()");
  if (static_cast<idx>(probXY.rows()) != X.size() || static_cast<idx>(probXY.cols()) != Y.size())
    throw exception::SizeMismatch("clara::cor()");
  return cov(probXY, X, Y) / (sigma(marginalX(probXY), X) * sigma(marginalX(probXY), Y));
}

}  // namespace clara

#endif  // !STATISTICS_H_
