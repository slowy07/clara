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

#ifndef CLASSFUNCTION_CODES_H_
#define CLASSFUNCTION_CODES_H_

#include "../functions.h"
#include "../internal/classFunction/singleton.h"
#include "../internal/util.h"
#include "exception.h"

namespace clara {

// codes class for handling different types of quantum error correction code
class Codes final : public internal::Singleton<const Codes> {
  friend class internal::Singleton<const Codes>;

 public:
  // enum to define different type of quantum error correction code
  enum class Type { FIVE_QUBIT, SEVEN_QUBIT_STEANE, NINE_QUBIT_SHOR };

 private:
  /**
   * default initialization of an object of const type ``const clara::Codes``,
   * requires a user-provided default constructor
   */
  Codes() = default;
  ~Codes() override = default;

 public:
  /**
   * @brief retrieve the codeword for the specified quantum error correction code
   *         type and index.
   * @param type the type of quantum error correction code
   * @param i the inde of the codeword to retrieve
   * @return the codeword as a ket (quantum state) for the specified type and index
   *
   * @throws NoCodeword exception if the specified type or index is invalid
   */
  ket codeword(Type type, idx i) const {
    ket result;
    switch (type) {
      case Type::FIVE_QUBIT:
        switch (i) {
          case 0:
            // define the codeword for the specified type and index
            result = (mket({0, 0, 0, 0, 0}) + mket({1, 0, 0, 1, 0}) + mket({0, 1, 0, 0, 1}) +
                      mket({1, 0, 1, 0, 0}) + mket({0, 1, 0, 1, 0}) - mket({1, 1, 0, 1, 1}) -
                      mket({0, 0, 1, 1, 0}) - mket({1, 1, 0, 0, 0}) - mket({1, 1, 1, 0, 1}) -
                      mket({0, 0, 0, 1, 1}) - mket({1, 1, 1, 1, 0}) - mket({0, 1, 1, 1, 1}) -
                      mket({1, 0, 0, 0, 1}) - mket({0, 1, 1, 0, 0}) - mket({1, 0, 1, 1, 1}) +
                      mket({0, 0, 1, 0, 1})) /
                     4.;
            break;
          case 1:
            // defined anoterh codeword for the specified type and index
            result = (mket({1, 1, 1, 1, 1}) + mket({0, 1, 1, 0, 1}) + mket({1, 0, 1, 1, 0}) +
                      mket({0, 1, 0, 1, 1}) + mket({1, 0, 1, 0, 1}) - mket({0, 0, 1, 0, 0}) -
                      mket({1, 1, 0, 0, 1}) - mket({0, 0, 1, 1, 1}) - mket({0, 0, 0, 1, 0}) -
                      mket({1, 1, 1, 0, 0}) - mket({0, 0, 0, 0, 1}) - mket({1, 0, 0, 0, 0}) -
                      mket({0, 1, 1, 1, 0}) - mket({1, 0, 0, 1, 1}) - mket({0, 1, 0, 0, 0}) +
                      mket({1, 1, 0, 1, 0})) /
                     4.;
            break;
          default:
            throw exception::NoCodeword("clara::Codes::codeword()");
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
            throw exception::NoCodeword("clara::Codes::codeword()");
        }
        break;
      case Type::NINE_QUBIT_SHOR:
        ket shora, shorb;
        // define quantum state for the shora and shorb kets
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
            // define qthe codeword using kronecker prodduct of shora ket
            result = kron(shora, kron(shora, shora)) / std::sqrt(8.);
            break;
          case 1:
            // defined the codeword using the kronecker product of shorb ket
            result = kron(shorb, kron(shorb, shorb)) / std::sqrt(8.);
            break;
          default:
            throw exception::NoCodeword("clara::Codes::codeword()");
        }
    }
    return result;
  }
};
}  // namespace clara

#endif
