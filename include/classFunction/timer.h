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

#ifndef CLASSFUNCTION_TIMER_H_
#define CLASSFUNCTION_TIMER_H_

#include <chrono>
#include <ostream>

#include "interface_display.h"
namespace clara {

/**
 * @class Timer
 * @brief class for measuring time intervals using the C++ chrono library
 *
 * the timer class provides functionally for measuring time intervals using
 * the c++ chrono library it can be used to measure the time taken for certain
 * operations or to benchmark code performance.
 *
 * @tparam T the type of duration to be used to measuring time
 * @param CLOCK_T clock type to be used
 */
template <typename T = std::chrono::duration<double>, typename CLOCK_T = std::chrono::steady_clock>
class Timer : public InterfaceDisplay {
 protected:
  typename CLOCK_T::time_point start_, end_;

 public:
  /**
   * @brief Construct an instance with the current time as starting point
   */
  Timer() noexcept : start_{CLOCK_T::now()}, end_{start_} {}

  ~Timer() override = default;

  /**
   * @brief Reset the chronometer
   * Reset the starting/ending point to the current time
   */
  void tic() noexcept {
    start_ = end_ = CLOCK_T::now();
    return *this;
  }

  /**
   * @brief Stops the chronometer
   * Set the current time as the ending point
   * @return Current instance
   */
  const Timer& toc() noexcept {
    end_ = CLOCK_T::now();
    return *this;
  }

  /**
   * @brief return the elapsed time as a duration of type U
   *
   * this function allow the caller to get the measured duration in any desired
   * resultion such as nanosecond, microsecond, millisecond, or second
   *
   * @tparam U type of duration to return (default = same as internal type T)
   * @return elapsed time as a duration of type U
   */
  template <typename U = T>
  U get_duration() const noexcept {
    return std::chrono::duration_cast<U>(end_ - start_);
  }

  /**
   * @return the elapsed time as floating-point count of internal duration unit
   *
   * useful when you want a simple numeric value instead of a duration object
   *
   * @return elapsed time in units of the internal duration type T as double
   */
  double tics() const noexcept {
    return static_cast<double>(std::chrono::duration_cast<T>(end_ - start_).count());
  }

  /**
   * @brief Get the duration in seconds
   * @return Number of seconds that have elapsed since the last call to tic() or reset()
   */
  double get_seconds() const noexcept { return get_duration().count(); }

 private:
  /**
   * @brief Override of clara::InterfaceDisplay::display()
   * Writes the output stream the number of seconds that have elapsed since the last call to
   * clara::Timer::toc()
   * @return The output stream
   */
  std::ostream& display(std::ostream& os) const override { return os << get_seconds(); }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_TIMER_H_
