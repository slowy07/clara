#ifndef CLASSFUNCTION_TIMER_H_
#define CLASSFUNCTION_TIMER_H_

#include <chrono>
#include <ostream>

#include "idisplay.h"
namespace clara {

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
   * @brief Duration specified by U
   * @return duration that passed between the reset and invocation
   * of clara::Timer::toc()
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
   * @brief clara::IDisplay::display() override
   * @return write to the output stream the number of tics
   * that passed between the reset and invocation of
   * clara::Timer::toc()
   */
  std::ostream& display(std::ostream& os) const override { return os << tics(); }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_TIMER_H_
