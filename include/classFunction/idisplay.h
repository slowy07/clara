#ifndef CLASSFUNCTION_IDISPLAY_H_
#define CLASSFUNCTION_IDISPLAY_H_

#include <ostream>
class IDisplay {
 private:
  virtual std::ostream& display(std::ostream& os) const = 0;

 public:
  IDisplay() = default;
  IDisplay(const IDisplay&) = default;
  IDisplay(IDisplay&&) = default;
  IDisplay& operator=(const IDisplay&) = default;
  virtual ~IDisplay() = default;
  friend inline std::ostream& operator<<(std::ostream& os, const IDisplay& rhs) {
    return rhs.display(os);
  }
};

#endif  // !CLASSFUNCTION_IDISPLAY_H_
