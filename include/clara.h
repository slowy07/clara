#ifndef CLARA_H_
#define CLARA_H_

#if (__GNUC__ && !__CLANG__)
#define CLARA_UNUSED_ __attribute__((unused))
#else
#define CLARA_UNUSED_
#endif  // (__GNUC__ && !__CLANG__)

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

// preprocessor macro
#include "macros.h"

// inter dependicies
#include "macros.h"

#include "types.h"
#include "classFunction/exception.h"
#include "constants.h"
#include "traits.h"
#include "classFunction/idisplay.h"
#include "internal/util.h"
#include "internal/classFunction/iomanip.h"
#include "input_output.h"

#include "internal/classFunction/singleton.h"
#include "classFunction/init.h"
#include "functions.h"
#include "classFunction/codes.h"
#include "classFunction/gates.h"
#include "classFunction/states.h"
#include "classFunction/random_devices.h"

#include "statistics.h"
#include "operations.h"
#include "entropies.h"
#include "entanglement.h"

#include "random.h"
#include "classFunction/timer.h"
#include "instruments.h"
#include "number_theory.h"

// order lib

namespace clara {
/**
 * @brief clara::Init const singleton
 * additional initialization/cleanups
 */
static const Init& init CLARA_UNUSED_ = Init::get_instance();

/**
 * @brief clara::Codes Signleton
 * initializes te code
 */
static const Codes& codes CLARA_UNUSED_ = Codes::get_instance();

/**
 * @brief clara::Gates const singleton
 * initializes the gates
 */
static const Gates& gt CLARA_UNUSED_ = Gates::get_instance();

/**
 * @brief clara::States const singleton
 * intializes the states
 */
static const States& st CLARA_UNUSED_ = States::get_instance();

/**
 * @brief Singleton
 * intialize the random device
 */
#ifdef NO_THREAD_LOCAL_
static RandomDevices& rdevs CLARA_UNUSED_ = RandomDevices::get_instance();
#else
thread_local static RandomDevices& rdevs CLARA_UNUSED_ = RandomDevices::get_thread_local_instance();
#endif  // DEBUG

}  // namespace clara

#endif  // !CLARA_H_
