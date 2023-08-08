#ifndef CLASSFUNCTION_IDISPLAY_H_
#define CLASSFUNCTION_IDISPLAY_H_

#include <ostream>
namespace clara {
/**
 * @class IDisplay
 * @brief Interface for displaying objects
 *
 * the IDisplay class provides as interface for object that need to be displayed in custom way
 * derived classes must override the `display()` function to perform the actual stream extraction
 * processing
 */
class IDisplay {
 private:
  /**
   * @brief pure virtual function for displaying objects
   *
   * this functions must be implemented in all derived classes to perform the actual
   * stream extraction processing. the function is automatically called by friend
   * inline function `operator<<`
   *
   * @param os the output stream to which the object should be displayed
   * @return the output stream after displaying the object
   */
  virtual std::ostream& display(std::ostream& os) const = 0;

 public:
  // default constructor
  IDisplay() = default;
  // copy constructor
  IDisplay(const IDisplay&) = default;
  // move constructor
  IDisplay(IDisplay&&) = default;
  // copy assignment operator
  IDisplay& operator=(const IDisplay&) = default;
  IDisplay& operator=(IDisplay&&) = default;
  // virtual destructor
  virtual ~IDisplay() = default;

  /**
   * @brief overloads the extraction operator to display objects
   *
   * @param os the output stream to which the object sould be displayed
   * @param rhs the object to be displayed
   * @return the output stream after displaying the object
   */
  friend inline std::ostream& operator<<(std::ostream& os, const IDisplay& rhs) {
    return rhs.display(os);
  }
};
}  // namespace clara

#endif  // !CLASSFUNCTION_IDISPLAY_H_
