include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_executable(
  clara_testing
  classFunction/codes_testing.cpp
  classFunction/gates_testing.cpp
  classFunction/idisplay_testing.cpp
  classFunction/random_devices_testing.cpp
  classFunction/reversible_testing.cpp
  classFunction/states_testing.cpp
  classFunction/timer_testing.cpp
  internal/classFunction/iomanip_testing.cpp
  internal/classFunction/singleton_testing.cpp
  testing_main.cpp
  traits.cpp
  input_output_testing.cpp
  operations.cpp
  random_clara.cpp)
target_link_libraries(clara_testing gmock)
