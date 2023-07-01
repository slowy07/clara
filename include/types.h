#ifndef TYPES_H_
#define TYPES_H_

#include <complex>
#include <eigen3/Eigen/Dense>

namespace clara {
namespace types {
typedef std::complex<double> cplx;

typedef Eigen::Matrix2cd cmat2;
typedef Eigen::Matrix3cd cmat3;
typedef Eigen::Matrix4cd cmat4;
typedef Eigen::MatrixXcd cmat;

typedef Eigen::Vector2cd ket2;
typedef Eigen::Vector3cd ket3;
typedef Eigen::Vector4cd ket4;
typedef Eigen::VectorXcd ket;

typedef Eigen::RowVector2cd bra2;
typedef Eigen::RowVector3cd bra3;
typedef Eigen::RowVector4cd bra4;
typedef Eigen::RowVectorXcd bra;

typedef Eigen::VectorXcd cvect;
typedef Eigen::VectorXi ivect;

}
}

#endif
