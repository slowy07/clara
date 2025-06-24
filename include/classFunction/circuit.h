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


#ifndef CLASSFUNCTION_CIRCUIT_H
#define CLASSFUNCTION_CIRCUIT_H

#include <vector>

#include "codes.h"
#include "gates.h"
#include "interface_display.h"

namespace clara {

class ClaraCircuitDescription : public InterfaceDisplay {
  const idx nq_;
  const idx nc_;
  const idx d_;
  std::vector<idx> measurement_steps_{};
  std::vector<bool> measured_;
  idx steps_cnt_;

 public:
  enum class GateType {
    NONE,
    SINGLE,
    TWO,
    THREE,
    CUSTOM,
    FAN,
    QFT,
    TFQ,
    SINGLE_CONTROL_SINGLE_TARGET,
    SINGLE_CONTROL_MULTIPLE_TARGET,
    MULTIPLE_CONTROL_SINGLE_TARGET,
    MULTIPLE_CONTROL_MULTIPLE_TARGET,
    CUSTOM_CONTROL,
    SINGLE_CCONTROL_SINGLE_TARGET,
    SINGLE_CCONTROL_MULTIPLE_TARGET,
    MULTIPLE_CCONTROL_MULTIPLE_TARGET,
    CUSTOM_CCONTROL,
  };

  enum class MeasureType {
    NONE,
    MEASURE_Z,
    MASURE_V_MANY,
  };

 public:
  class iterator {
    friend ClaraCircuitDescription;
    friend class ClaraCricuit;

    const ClaraCircuitDescription* qcd_{nullptr};

    struct value_type_ : public InterfaceDisplay {};
  };
};

}  // namespace clara

#endif  // !CLASSFUNCTION_CIRCUIT_H
