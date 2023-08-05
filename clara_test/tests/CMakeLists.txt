include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_executable(
  clara_testing
  classFunction/gates_testing.cpp
  classFunction/states_testing.cpp
  classFunction/timer_testing.cpp
  classFunction/codes_testing.cpp
  classFunction/random_devices_testing.cpp
  testing_main.cpp
  traits.cpp
  operations.cpp
  random_clara.cpp)
target_link_libraries(clara_testing gmock)
