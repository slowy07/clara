#ifndef TYPES_H_
#define TYPES_H_

#include <complex>
#include <eigen3/Eigen/Dense>

namespace clara {
namespace types {
typedef std::complex<double> cplx;

// complex matrix
typedef Eigen::MatrixXcd cmat;

// double matrix
typedef Eigen::MatrixXd dmat;

// integer Matrix
typedef Eigen::MatrixXi imat;

}
}

#endif
