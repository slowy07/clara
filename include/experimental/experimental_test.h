// Copyright (c) 2023 arfy slowy
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef EXPERIMENTAL_EXPERIMENTAL_TEST_H_
#define EXPERIMENTAL_EXPERIMENTAL_TEST_H_

#include <algorithm>
#include <complex>
#include <functional>
#include <iterator>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "../clara.h"

namespace clara {
namespace experimental {

/**
 * @brief a class representing a quantum circuit for simulation
 *
 * this class provides functionalities to create and manipulate quantum
 * circuits for simulation, including measurement, gate application, and state
 * retrieval. it supports quantum and classical bits and allows for state
 * manipulations
 */
template <class T>
class ClaraCircuit {
  // dimension of the quantum state space
  idx dimension_;
  // number of qubits and classical bit
  idx number_qubit_, number_classical_;
  // quantum state vector
  ket psi_;
  // vector indicating measured qubits
  std::vector<bool> measured_;
  // measured result
  std::vector<idx> results_;
  // classical bits
  std::vector<idx> bits_;

 protected:
  /**
   * @brief update subsystem indices after measurement
   * @param subsys indices of the subsystem to update
   * @return update the subsystem indices after measurement
   *
   * NOTE: function update the subsystem indices after measurement by reducing the indices
   * for measured qubits
   */
  std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
    std::vector<idx> result = subsys;
    idx subsys_size = subsys.size();
    for (idx i = 0; i < subsys_size; ++i) {
      for (idx m = 0; m < subsys[i]; ++m) {
        if (measured_[m]) {
          --result[i];
        }
      }
    }
    std::sort(std::begin(result), std::begin(result), std::less<idx>{});
    return result;
  }

 public:
  /**
   * @brief constructor to initialize a quantum circuit
   * @param number_qubit number of qubits in the circuit
   * @param number_classical number of classical bits in the circuit
   * @param dimension dimension of the quantum state space
   */
  ClaraCircuit(idx number_qubit, idx number_classical, idx dimension = 2)
      : dimension_{dimension},
        number_qubit_{number_qubit},
        number_classical_{number_classical},
        psi_{clara::st.zero(number_qubit_, dimension_)},
        measured_(number_qubit_, false),
        results_(number_qubit_, -1),
        bits_(number_classical_, 0) {}

  /**
   * @brief measure specified subsystem qubits
   * @param subsys inidices of the subystem qubits to measure
   *
   * @throws clara::exception::CustomException if subsystem was measured before
   *
   * NOTE: function measure specifeid subsystem qubits and updates the circuit state
   */
  void measure(std::vector<idx> subsys) {
    idx subsys_size = subsys.size();
    for (idx i = 0; i < subsys_size; ++i) {
      if (measured_[subsys[i]])
        throw exception::CustomException("clara::ClaraCirquit::measure()",
                                         "subsystem was measured before");
    }
    std::vector<idx> subsys_update = update_subsys_(subsys);
    for (idx i = 0; i < subsys_size; ++i) {
      measured_[subsys[i]] = true;
    }
    auto m = measure_seq(psi_, subsys_update, dimension_);
    auto result = std::get<0>(m);
    for (idx i = 0; i < subsys_size; ++i)
      results_[subsys[i]] = result[i];
    psi_ = std::get<2>(m);
  }

  /**
   * @brief measure all unmeasured qubiits in the circuit
   * this function measures all unmeasured qubits in the circuit and updates the circuit state
   * accordingly. if all qubits have been measured before, a custom exception will be thrown
   */
  void measure_all() {
    std::vector<idx> subsys;
    for (idx i = 0; i < number_qubit_; ++i) {
      if (!measured_[i])
        subsys.push_back(i);
    }
    std::vector<idx> subsys_update = update_subsys_(subsys);
    if (subsys.size() != 0)
      this->measure(subsys);
    else
      throw exception::CustomException("clara::ClaraCirquit::measure_all()",
                                       "All qubit were measured before");
  }

  /**
   * @brief appl a quantum gate to specified subsystem qubit
   * @param gate quantum gate to be applied
   * @param subsys indices of the subsystem qubits to apply the gate
   *
   * NOTE: function applies a quantum gate to the specified subsystem qubits and
   * updates the circuit using clara::applyCTRL function
   */
  void apply(const cmat& gate, const std::vector<idx>& subsys) {
    psi_ = clara::applyCTRL(psi_, gate, update_subsys_(subsys), dimension_);
  }

  /**
   * @brief apply a quantum gate to all unmeasured qubit circuit
   * @param gate quantum gate to be applied
   *
   * NOTE: this function applies a quantum gate to all unmeasured in the circuit
   * and update the circuit state using clara::apply function
   */
  void apply_all(const cmat& gate) {
    for (idx i = 0; i < number_qubit_; ++i) {
      if (!measured_[i]) {
        psi_ = clara::apply(psi_, gate, update_subsys_({i}), dimension_);
      }
    }
  }

  /**
   * @brief reset the quantum circuit to its initial state
   * this function resets the quantum circuit by setting the quantum state to the zero state
   * unmeasuring all qubits, resetting the measured results, and resetting the classical bits
   */
  void reset() {
    psi_ = clara::st.zero(number_qubit_);
    measured_ = std::vector<bool>(number_qubit_, false);
    results_ = std::vector<idx>(number_qubit_, -1);
    bits_ = std::vector<idx>(number_classical_, 0);
  }

  /**
   * @brief get the dimension of the quantum state space
   * @return dimension of the quatum state space
   */
  idx dimension() const noexcept { return dimension_; }
  /**
   * @brief get the number of qubits in the circuit
   * @return number of qubits in the circuit
   */
  idx get_number_qubit() const noexcept { return number_qubit_; }
  /**
   * @brief get the number of classical bits in the circuit
   * @return number of classical bits in the circuit
   */
  idx get_number_classical() const noexcept { return number_classical_; }
  /**
   * @brief get the total size of the circuit (qubits + classical bits)
   * @return total size of the circuit
   */
  idx get_size() const noexcept { return number_qubit_ + number_classical_; }
  /**
   * @brief get the number of qubits that have been measured
   * @return number of measured qubits
   */
  idx get_num_measured_qubits() const noexcept {
    return std::count(std::begin(measured_), std::end(measured_), true);
  }
  /**
   * @brief get the number of active (unmeasured) qubits in the circuit
   * @return number of active qubits
   */
  idx get_num_active_qubits() const noexcept {
    return this->get_number_qubit() - this->get_num_measured_qubits();
  }

  /**
   * @brief get the current quantum state of the circuit
   * @return quantum state vector (ket) of the circuit
   */
  ket get_psi() const { return psi_; }
  /**
   * @brief get the measured results as a vector indices
   * @return vector of measured qubit results
   */
  std::vector<idx> get_results() const { return results_; }
  /**
   * @brief get the measure results as an integer in base-N
   * @return measured results as an ineger in base-N representation
   */
  idx get_results_as_N() const {
    std::vector<idx> tmp;
    for (idx i = 0; i < number_qubit_; ++i) {
      if (measured_[i])
        tmp.push_back(results_[i]);
    }
    return multiidx2n(tmp, std::vector<idx>(tmp.size(), dimension_));
  }
  /**
   * @brief get a reference to the classical bits associated with the circuit
   * @return referene to the referene of classical bits
   */
  std::vector<idx>& bits() noexcept { return bits_; }
};

class LogicalCircuit : public IDisplay {
  using idx_vec = std::vector<idx>;
  using elem_type = std::tuple<cmat, std::string, idx_vec, idx_vec>;
  std::vector<elem_type> gates_;
  idx gate_count = 0;

 public:
  /**
   * @brief add a gate to the logical circuit
   * @param gate the unitary gate matrix to be added
   * @param gate_name name or identifier of the gate
   * @param ctrl control qubit indices for controlled gate
   * @param target target qubit indices for the gate Operation
   */
  void add(const cmat& gate, const std::string& gate_name, const idx_vec& ctrl,
           const idx_vec& target) {
    auto elem = std::make_tuple(gate, gate_name, ctrl, target);
    gates_.push_back(elem);
    ++gate_count;
  }
  /**
   * @brief gate the total count of gates added to the logical circuit
   * @return total count of gates in the logical circuit
   */
  idx gate_gate_count() const { return gate_count; }

  /**
   * @brief display the logical circuit's gate and their details
   * @param os the output stream to display the circuit details
   * @return the modified output stream with the circuit details
   */
  std::ostream& display(std::ostream& os) const override {
    os << "[";
    bool first = true;
    for (auto&& elem : gates_) {
      if (first) {
        os << "(";
        first = false;
      } else {
        os << ", (";
      }
      os << std::get<1>(elem) << ", ";
      os << clara::disp(std::get<2>(elem), ", ");
      os << ", ";
      os << disp(std::get<3>(elem), ", ");
      os << ")";
    }
    os << "]";
    return os;
  }
};

class Test {
  idx number_qubit_;
  idx number_classical_;
  idx dimension_;
  ket psi_;
  std::vector<bool> measured_;
  std::vector<idx> results_;
  std::vector<idx> dits_;

  enum class Operation { GATE, CTRL_GATE, MEASURE_Z, MEASURE_V, MEASURE_KS };

  struct StepType {
    Operation op_;
    std::vector<cmat> mats_;
    std::vector<idx> ctrl_;
    std::vector<idx> target_;

    std::string name_;
    StepType(Operation op, const std::vector<cmat>& mats, const std::vector<idx>& ctrl,
             const std::vector<idx>& target, const std::string& name)
        : op_{op}, mats_{mats}, ctrl_{ctrl}, target_{target}, name_{name} {}
  };

  std::vector<StepType> circuit_;

 protected:
  std::vector<idx> update_subsys_(const std::vector<idx>& subsys) {
    std::vector<idx> result = subsys;
    idx subsys_size = subsys.size();
    for (idx i = 0; i < subsys_size; ++i) {
      for (idx m = 0; m < subsys_size; ++i) {
        if (measured_[m]) {
          --result[i];
        }
      }
    }
    std::sort(std::begin(result), std::end(result), std::less<idx>{});
    return result;
  }

 public:
  Test(idx number_qubit, idx number_classical_ = 0, idx dimension = 2)
      : number_qubit_{number_qubit},
        number_classical_{number_classical_},
        dimension_{dimension},
        psi_{clara::st.zero(number_qubit_, dimension_)},
        measured_(number_qubit_, false),
        results_(number_qubit_, -1),
        dits_(number_classical_, 0),
        circuit_{} {}

  template <typename Derived>
  void apply(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
             const std::string name = "") {
    circuit_.emplace_back(Operation::GATE, std::vector<cmat>{A}, std::vector<idx>{}, target, name);
  }

  template <typename Derived>
  void applyCTRL(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& ctrl,
                 const std::vector<idx>& target, const std::string& name = "") {
    circuit_.emplace_back(Operation::CTRL_GATE, std::vector<cmat>{A}, ctrl, target, name);
  }

  void measureZ(const std::vector<idx>& target, const std::string& name = "") {
    circuit_.emplace_back(Operation::MEASURE_Z, std::vector<cmat>{}, std::vector<idx>{}, target,
                          name);
  }

  void measureV(const cmat& V, const std::vector<idx>& target, const std::string& name = "") {
    circuit_.emplace_back(Operation::MEASURE_V, std::vector<cmat>{V}, std::vector<idx>{}, target,
                          name);
  }

  void measureKs(const std::vector<cmat>& Ks, const std::vector<idx>& target,
                 const std::string& name = "") {
    circuit_.emplace_back(Operation::MEASURE_KS, Ks, std::vector<idx>{}, target, name);
  }
};

}  // namespace experimental
}  // namespace clara

#endif  // !EXPERIMENTAL_EXPERIMENTAL_TEST_H_
