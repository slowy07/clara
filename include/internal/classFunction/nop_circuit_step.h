#ifndef INTERNAL_CLASSFUNCTION_NOP_CIRCUIT_STEP_H_
#define INTERNAL_CLASSFUNCTION_NOP_CIRCUIT_STEP_H_

#include "../../classFunction/interface_display.h"

namespace clara {
namespace internal {

/**
 * @brief representing "no operation (NOP)" step in circuit simulation
 *
 * this class model step that does nothing when get execute, it useful as placeholder,
 * default step, or representing idle state in sequence operation
 */
struct CircuitNOPStep : InterfaceDisplay {
  /**
   * @brief compare two CircuitNOPStep object for equality
   *
   * since NOP step have no internal state, all instance are considered are equals
   *
   * @param other the other CircuitNOPStep to compare
   * @return always return true
   */
  bool operator==(const CircuitNOPStep&) const noexcept { return true; }

  /**
   * @brief compare two CircuitNOPStep object for inequality
   *
   * relies that the == operator to determine if the object are not equal
   *
   * @param rhs right-hand side CircuitNOPStep to compare with
   * @return false because all NOP step are considered equals
   */
  bool operator!=(const CircuitNOPStep& rhs) const noexcept { return !(*this == rhs); }

 private:
  /**
   * @brief output string representation of this step to a streams
   *
   * implements the pure virtual display() method from the InterfaceDisplay base class
   * this will allow CircuitNOPStep to be printed or logged consistenly with other step
   *
   * @param os the output stream to write to
   * @return a reference to the modified output stream for chaining
   */
  std::ostream& display(std::ostream& os) const override {
    os << "NOP";  // writing NOP on stream
    return os;
  }
};
}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_NOP_CIRCUIT_STEP_H_
