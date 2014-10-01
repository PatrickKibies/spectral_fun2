cmake_minimum_required(VERSION 2.8)
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR/cmake}")
project(spectral_fun2)

SET(FFTW_INCLUDES "/home/patrick/lib/fftw3/include")
SET(FFTW_LIBRARIES "/home/patrick/lib/fftw3/lib")

#find_package(boost)

add_executable(spectral_fun2 main.cpp)

target_link_libraries(spectral_fun2 fftw3)
target_link_libraries(spectral_fun2 m)

install(TARGETS spectral_fun2 RUNTIME DESTINATION bin)
