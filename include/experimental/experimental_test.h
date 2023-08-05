#ifndef EXPERIMENTAL_EXPERIMENTAL_TEST_H_
#define EXPERIMENTAL_EXPERIMENTAL_TEST_H_

#include <algorithm>
#include <complex>
#include <functional>
#include <iterator>
#include <ostream>
#include <string>
#include <vector>

#include "../classFunction/gates.h"
#include "../classFunction/idisplay.h"

namespace clara {
namespace experimental {

class ClaraCircuit : public IDisplay {
  idx number_qudits_;
  idx number_classical_bits_;
  idx dimension_;
  ket psi_;
  std::vector<idx> dits_;
  std::vector<idx> measurement_steps_;

  enum class GateType {
    SINGLE,
    TWO,
    THREE,
    FAN,
    CUSTOM,
    SINGLE_CTRL_SINGLE_TARGET,
    SINGLE_CTRL_MULTIPLE_TARGET,
    MULTIPLE_CTRL_SINGLE_TARGET,
    MULTIPLE_CTRL_MULTIPLE_TARGET,
    CUSTOM_CTRL,
    SINGLE_cCTRL_SINGLE_TARGET,
    SINGLE_cCTRL_MULTIPLE_TARGET,
    MULTIPLE_cCTRL_SINGLE_TARGET,
    MULTIPLE_cCTRL_MULTIPLE_TARGET,
    CUSTOM_cCTRL,
    NONE
  };

  enum class MeasureType { MEASURE_Z, MEASURE_V, MEASURE_KS, NONE };

  struct GateStep {
    GateType gate_type_ = GateType::NONE;
    cmat gate_;
    std::vector<idx> ctrl_;
    std::vector<idx> target_;
    std::string name_;

    GateStep() = default;
    GateStep(GateType gate_type, const cmat& gate, const std::vector<idx>& ctrl,
             const std::vector<idx>& target, const std::string& name = "")
        : gate_type_{gate_type}, gate_{gate}, ctrl_{ctrl}, target_{target}, name_{name} {}
  };

  struct MeasureStep {
    MeasureType measurement_type_ = MeasureType::NONE;
    std::vector<cmat> mats_;
    std::vector<idx> target_;
    std::string name_;

    MeasureStep() = default;
    MeasureStep(MeasureType measurement_type, const std::vector<cmat>& mats,
                const std::vector<idx>& target, const std::string& name = "")
        : measurement_type_{measurement_type}, mats_{mats}, target_{target}, name_{name} {}
  };

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
    std::sort(std::begin(result), std::end(result), std::less<idx>{});
    return result;
  }

 public:
  std::vector<bool> measured_;
  ClaraCircuit(idx number_qudit, idx number_classical_bits = 0, idx dimension = 2,
               const std::string& name = "");
  void apply(const cmat& U, idx i, const std::string& name = "");
  void measureZ(idx i, const std::string& name = "");
  std::ostream& display(std::ostream& os) const override;

  friend std::ostream& operator<<(std::ostream& os, const GateType& gate_type);
  friend std::ostream& operator<<(std::ostream& os, const MeasureType& measurement_type);
};

}  // namespace experimental
}  // namespace clara

#endif  // !EXPERIMENTAL_EXPERIMENTAL_TEST_H_
