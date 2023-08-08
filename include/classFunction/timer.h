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
template <typename T = std::chrono::duration<double>, typename CLOCK_T = std::chrono::steady_clock>
class Timer : public IDisplay {
 protected:
  typename CLOCK_T::time_point start_, end_;

 public:
  /**
   * @brief construct an instance with the current time
   * as starting point
   */
  Timer() noexcept : start_{CLOCK_T::now()}, end_{start_} {}

  /**
   * @brief reset the chronometer
   * reset the starting/ending point to the current time
   */
  void tic() noexcept { start_ = end_ = CLOCK_T::now(); }

  /**
   * @brief stops the chronometer
   * set the current time as the ending point
   * @return current instance
   */
  const Timer& toc() noexcept {
    end_ = CLOCK_T::now();
    return *this;
  }
  /**
   * @brief time passsed in the duration specified by T
   * @return number of tics (specified by T) that passed between the
   * instantiation/rest and invocation of clara::Timer::toc()
   */
  double tics() const noexcept { return std::chrono::duration_cast<T>(end_ - start_).count(); }
  /**
   * @brief get the duration specified by U
   *
   * get the duration that passed between the reset and invocation of clara::Timer::toc()
   * @tparam U the type of duration to be returned
   * @return Duration that passed between that reset and invocation clara::Timer::toc()
   */
  template <typename U = T>
  U get_duration() const noexcept {
    return std::chrono::duration_cast<U>(end_ - start_);
  }
  /**
   * @brief default copy constructor
   */
  Timer(const Timer&) = default;
  /**
   * @brief default move constructor
   */
  Timer(Timer&&) = default;
  /**
   * @brief default copy assignment operator
   */
  Timer& operator=(const Timer&) = default;
  /**
   * @brief default move assignment operator
   */
  Timer& operator=(Timer&&) = default;
  /**
   * @brief default virtual destructor
   */
  virtual ~Timer() = default;

 private:
  /**
   * @brief override of clara::IDisplay::display()
   * writes the output strem the number of tics that passed between the reset and invocation of
   * clara::TImer::toc()
   * @return the output stream
   */
  std::ostream& display(std::ostream& os) const override { return os << tics(); }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_TIMER_H_
