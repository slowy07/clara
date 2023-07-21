include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_executable(
  clara_testing classFunction/gates.cpp classFunction/states.cpp
                classFunction/timer.cpp testing_main.cpp traits.cpp)
target_link_libraries(clara_testing gmock)
