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
