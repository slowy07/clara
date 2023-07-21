#ifndef CLASSFUNCTION_STATES_H_
#define CLASSFUNCTION_STATES_H_

#include <cmath>
#include <vector>

#include "../internal/classFunction/singleton.h"
#include "codes.h"
#include "exception.h"

namespace clara {

/**
 * @class clara::States
 * @brief const singleton class that implements most commonly used states
 */

class States final : public internal::Singleton<const States> {
  friend class internal::Singleton<const States>;

 public:
  ket x0{ket::Zero(2)};
  ket x1{ket::Zero(2)};
  ket y0{ket::Zero(2)};
  ket y1{ket::Zero(2)};
  ket z0{ket::Zero(2)};
  ket z1{ket::Zero(2)};

  cmat px0{cmat::Zero(2, 2)};
  // < Projector onto the Pauli Sigma-X 0-eigenstate |+><+|
  cmat px1{cmat::Zero(2, 2)};
  // < Projector onto the Pauli Sigma-X 1-eigenstate |-><-|
  cmat py0{cmat::Zero(2, 2)};
  // < Projector onto the Pauli Sigma-Y 0-eigenstate |y+><y+|
  cmat py1{cmat::Zero(2, 2)};
  // < Projector onto the Pauli Sigma-Y 1-eigenstate |y-><y-|
  cmat pz0{cmat::Zero(2, 2)};
  // < Projector onto the Pauli Sigma-Z 0-eigenstate |0><0|
  cmat pz1{cmat::Zero(2, 2)};
  // Bell states
  ket b00{ket::Zero(4)};
  ket b01{ket::Zero(4)};
  ket b10{ket::Zero(4)};
  ket b11{ket::Zero(4)};

  cmat pb00{cmat::Zero(4, 4)};
  cmat pb01{cmat::Zero(4, 4)};
  cmat pb10{cmat::Zero(4, 4)};
  cmat pb11{cmat::Zero(4, 4)};
  // W and GHZ states
  ket GHZ{ket::Zero(8)};
  ket W{ket::Zero(8)};

  // projectors onto GHZ and W
  cmat pGHZ{cmat::Zero(8, 8)};
  cmat pW{cmat::Zero(8, 8)};

  /**
   * @brief maximally entangled state of 2 qudits
   *  \f$\frac{1}{\sqrt{d}}\sum_{j=0}^{d-1}|jj\rangle\f$ of 2 qudits
   */
  ket mes(idx d = 2) const {
    // check valid dims
    if (d == 0)
      throw exception::DimsInvalid("clara::States::mes()");
    ket psi = mket({0, 0}, {d, d});
    for (idx i = 1; i < d; ++i) {
      psi += mket({i, i}, {d, d});
    }
    return psi / std::sqrt(d);
  }

  /**
   * @brief zero state of n qudits
   * @return zero state \f$|0\rangle^{\otimes n}\f$ of n qudits
   */
  ket zero(idx n, idx d = 2) const {
    if (n == 0)
      throw exception::OutOfRange("clara::States::zero()");
    if (d == 0)
      throw exception::DimsInvalid("clara::States::zero()");
    idx D = static_cast<idx>(std::pow(d, n));
    ket result = ket::Zero(D);
    result(0) = 1;

    return result;
  }

  /**
   * @brief one state of n qudits
   * @return one statem f$|1\rangle^{\otimes n}\f$ of qudits
   */
  ket one(idx n, idx d = 2) const {
    if (n == 0)
      throw exception::OutOfRange("clara::States::one()");
    if (d == 0)
      throw exception::DimsInvalid("clara::States::one()");
    ket result = ket::Zero(static_cast<ket::Index>(std::pow(d, n)));
    result(multiidx2n(std::vector<idx>(n, 1), std::vector<idx>(n, d))) = 1;
    return result;
  }

  /**
   * @brief \f$|j\rangle^{\otimes n}\f$ state of n qudits
   * @return \f$|j\rangle^{\otimes n}\f$ state of n qudits
   */
  ket jn(idx j, idx n, idx d = 2) const {
    if (n == 0)
      throw exception::OutOfRange("clara::States::jn()");
    if (j >= d)
      throw exception::SubsysMismatchdims("clara::States::jn()");

    if (j >= d)
      throw exception::DimsInvalid("clara::States::jn()");

    ket result = ket::Zero(static_cast<ket::Index>(std::pow(d, n)));
    result(multiidx2n(std::vector<idx>(n, j), std::vector<idx>(n, d))) = 1;
    return result;
  }

  /**
  * @brief plus state of n qubits
  * @return plus state \f$|+\rangle^{\otimes n}\f$ of qubits
*/
  ket plus(idx n) const {
    if (n == 0)
      throw exception::OutOfRange("clara::States::plus()");
    idx D = static_cast<idx>(std::pow(2, n));
    ket result = ket::Ones(D);

    return result / std::sqrt(D);
  }

  /**
  * @brief minus state of n qubits
  * @return minus state \f$|-\rangle^{\otimes n}\f$ of n qubits
*/
  ket minus(idx n) const {
    if (n == 0)
      throw exception::OutOfRange("clara::States::minus()");
    return kronpow(this -> x1, n);
  }

 private:
  States() {
    x0 << 1 / std::sqrt(2.), 1 / std::sqrt(2.);
    x1 << 1 / std::sqrt(2.), -1 / std::sqrt(2.);
    y0 << 1 / std::sqrt(2.), 1_i / std::sqrt(2.);
    y1 << 1 / std::sqrt(2.), -1_i / std::sqrt(2.);
    z0 << 1, 0;
    z1 << 0, 1;
    px0 = x0 * x0.adjoint();
    px1 = x1 * x1.adjoint();
    py0 = y0 * y0.adjoint();
    py1 = y1 * y1.adjoint();
    pz0 = z0 * z0.adjoint();
    pz1 = z1 * z1.adjoint();

    b00 << 1 / std::sqrt(2.), 0, 0, 1 / std::sqrt(2.);
    // (|00> + |11>) / sqrt(2)
    b01 << 0, 1 / std::sqrt(2.), 1 / std::sqrt(2.), 0;
    // (|01> + |10>) / sqrt(2)
    b10 << 1 / std::sqrt(2.), 0, 0, -1 / std::sqrt(2.);
    // (|00> - |11>) / sqrt(2)
    b11 << 0, 1 / std::sqrt(2.), -1 / std::sqrt(2.), 0;
    // (|01> - |10>) / sqrt(2)

    pb00 = b00 * b00.adjoint();
    pb01 = b01 * b01.adjoint();
    pb10 = b10 * b10.adjoint();
    pb11 = b11 * b11.adjoint();

    GHZ << 1, 0, 0, 0, 0, 0, 0, 1;
    GHZ = GHZ / std::sqrt(2.);
    W << 0, 1, 1, 0, 1, 0, 0, 0;
    W = W / std::sqrt(3.);

    pGHZ = GHZ * GHZ.adjoint();
    pW = W * W.adjoint();
  }
  ~States() = default;
};

}  // namespace clara

#endif  // !CLASSFUNCTION_STATES_H_
