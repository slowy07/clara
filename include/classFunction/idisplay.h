#ifndef CLASSFUNCTION_IDISPLAY_H_
#define CLASSFUNCTION_IDISPLAY_H_

#include <ostream>
namespace clara {
class IDisplay {
 private:
  /**
   * @brief this function must be implemented in all derived classes
   * the derived classes overrided this function to perform the actual
   * stream extraction processing. the function is automatically called
   * by the friend inline function
   */
  virtual std::ostream& display(std::ostream& os) const = 0;

 public:
  IDisplay() = default;
  IDisplay(const IDisplay&) = default;
  IDisplay(IDisplay&&) = default;
  IDisplay& operator=(const IDisplay&) = default;
  virtual ~IDisplay() = default;

  /**
   * @brief overloads the extraction operator
   */
  friend inline std::ostream& operator<<(std::ostream& os, const IDisplay& rhs) {
    return rhs.display(os);
  }
};
}  // namespace clara

#endif  // !CLASSFUNCTION_IDISPLAY_H_
