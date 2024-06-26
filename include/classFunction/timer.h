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

#include "idisplay.h"
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
template <typename Clock = std::chrono::steady_clock, typename Duration = std::chrono::duration<double>>
class Timer : public IDisplay {
 protected:
  typename Clock::time_point start_, end_;

 public:
  /**
   * @brief Construct an instance with the current time as starting point
   */
  Timer() noexcept : start_{Clock::now()}, end_{start_} {}

  /**
   * @brief Reset the chronometer
   * Reset the starting/ending point to the current time
   */
  void tic() noexcept { start_ = end_ = Clock::now(); }

  /**
   * @brief Stops the chronometer
   * Set the current time as the ending point
   * @return Current instance
   */
  const Timer& toc() noexcept {
    end_ = Clock::now();
    return *this;
  }

  /**
   * @brief Get the duration of the time interval
   * @return Duration that has elapsed since the last call to tic() or reset()
   */
  Duration get_duration() const noexcept {
    return std::chrono::duration_cast<Duration>(end_ - start_);
  }

  /**
   * @brief Get the duration in seconds
   * @return Number of seconds that have elapsed since the last call to tic() or reset()
   */
  double get_seconds() const noexcept { return get_duration().count(); }

  /**
   * @brief Get the duration in milliseconds
   * @return Number of milliseconds that have elapsed since the last call to tic() or reset()
   */
  double get_milliseconds() const noexcept {
    return std::chrono::duration_cast<std::chrono::milliseconds>(get_duration()).count();
  }

  /**
   * @brief Default copy constructor
   */
  Timer(const Timer&) = default;

  /**
   * @brief Default move constructor
   */
  Timer(Timer&&) = default;

  /**
   * @brief Default copy assignment operator
   */
  Timer& operator=(const Timer&) = default;

  /**
   * @brief Default move assignment operator
   */
  Timer& operator=(Timer&&) = default;

  /**
   * @brief Default virtual destructor
   */
  virtual ~Timer() = default;

 private:
  /**
   * @brief Override of clara::IDisplay::display()
   * Writes the output stream the number of seconds that have elapsed since the last call to
   * clara::Timer::toc()
   * @return The output stream
   */
  std::ostream& display(std::ostream& os) const override { return os << get_seconds(); }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_TIMER_H_
