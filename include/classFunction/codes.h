#ifndef CLASSFUNCTION_CODES_H_
#define CLASSFUNCTION_CODES_H_

#include "../functions.h"
#include "../internal/classFunction/singleton.h"
#include "../internal/util.h"
#include "exception.h"

namespace clara {
class Codes final : public internal::Singleton<const Codes> {
  friend class internal::Singleton<const Codes>;

 public:
  enum class Type { FIVE_QUBIT = 1, SEVEN_QUBIT_STEANE, NINE_QUBIT_SHOR };

 private:
  /**
   * default initialization of an object of const type ``const clara::Codes``,
   * requires a user-provided default constructor
   */
  Codes() {}
  ~Codes() = default;

 public:
  ket codeword(Type type, idx i) const {
    ket result;
    switch (type) {
      case Type::FIVE_QUBIT:
        switch (i) {
          case 0:
            result = (mket({0, 0, 0, 0, 0}) + mket({1, 0, 0, 1, 0}) + mket({0, 1, 0, 0, 1}) +
                      mket({1, 0, 1, 0, 0}) + mket({0, 1, 0, 1, 0}) - mket({1, 1, 0, 1, 1}) -
                      mket({0, 0, 1, 1, 0}) - mket({1, 1, 0, 0, 0}) - mket({1, 1, 1, 0, 1}) -
                      mket({0, 0, 0, 1, 1}) - mket({1, 1, 1, 1, 0}) - mket({0, 1, 1, 1, 1}) -
                      mket({1, 0, 0, 0, 1}) - mket({0, 1, 1, 0, 0}) - mket({1, 0, 1, 1, 1}) +
                      mket({0, 0, 1, 0, 1})) /
                     4.;
            break;
          case 1:
            result = (mket({1, 1, 1, 1, 1}) + mket({0, 1, 1, 0, 1}) + mket({1, 0, 1, 1, 0}) +
                      mket({0, 1, 0, 1, 1}) + mket({1, 0, 1, 0, 1}) - mket({0, 0, 1, 0, 0}) -
                      mket({1, 1, 0, 0, 1}) - mket({0, 0, 1, 1, 1}) - mket({0, 0, 0, 1, 0}) -
                      mket({1, 1, 1, 0, 0}) - mket({0, 0, 0, 0, 1}) - mket({1, 0, 0, 0, 0}) -
                      mket({0, 1, 1, 1, 0}) - mket({1, 0, 0, 1, 1}) - mket({0, 1, 0, 0, 0}) +
                      mket({1, 1, 0, 1, 0})) /
                     4.;
            break;
          default:
            throw Exception("clara::Codes::codeword()", Exception::Type::NO_CODEWORD);
        }
        break;
      case Type::SEVEN_QUBIT_STEANE:
        switch (i) {
          case 0:
            result = (mket({0, 0, 0, 0, 0, 0, 0}) + mket({1, 0, 1, 0, 1, 0, 1}) +
                      mket({0, 1, 1, 0, 0, 1, 1}) + mket({1, 1, 0, 0, 1, 1, 0}) +
                      mket({0, 0, 0, 1, 1, 1, 1}) + mket({1, 0, 1, 1, 0, 1, 0}) +
                      mket({0, 1, 1, 1, 1, 0, 0}) + mket({1, 1, 0, 1, 0, 0, 1})) /
                     std::sqrt(8.);
            break;
          case 1:
            result = (mket({1, 1, 1, 1, 1, 1, 1}) + mket({0, 1, 0, 1, 0, 1, 0}) +
                      mket({1, 0, 0, 1, 1, 0, 0}) + mket({0, 0, 1, 1, 0, 0, 1}) +
                      mket({1, 1, 1, 0, 0, 0, 0}) + mket({0, 1, 0, 0, 1, 0, 1}) +
                      mket({1, 0, 0, 0, 0, 1, 1}) + mket({0, 0, 1, 0, 1, 1, 0})) /
                     std::sqrt(8.);
            break;
          default:
            throw Exception("clara::Codes::codeword()", Exception::Type::NO_CODEWORD);
        }
        break;
      case Type::NINE_QUBIT_SHOR:
        ket shora, shorb;
        shora = mket({0, 0, 0}) + mket({
                                      1,
                                      1,
                                      1,
                                  });
        shorb = mket({0, 0, 0}) - mket({
                                      1,
                                      1,
                                      1,
                                  });
        switch (i) {
          case 0:
            result = kron(shora, kron(shora, shora)) / std::sqrt(8.);
            break;
          case 1:
            result = kron(shorb, kron(shorb, shorb)) / std::sqrt(8.);
            break;
          default:
            throw Exception("clara::Codes::codeword()", Exception::Type::NO_CODEWORD);
        }
    }
    return result;
  }
};
}  // namespace clara

#endif
