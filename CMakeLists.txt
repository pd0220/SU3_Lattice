cmake_minimum_required(VERSION 3.0.0)
project(SU3_Lattice VERSION 0.1.0)

include(CTest)
enable_testing()

# GSL
# INCLUDE_DIRECTORIES("C:\\DEV\\vcpkg\\installed\\x86-windows\\include")
# INCLUDE_DIRECTORIES("C:\\DEV\\vcpkg\\installed\\x86-windows\\lib")
# EIGEN3
SET( EIGEN3_INCLUDE_DIR "C:\\Program\ Files\\eigen-3.4.0")
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

if (MSVC)
  string(REGEX REPLACE "/W[0-9]" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif (MSVC)

add_executable(SU3_Lattice main.cpp)

set_target_properties(SU3_Lattice PROPERTIES CXX_STANDARD 17
                                CXX_STANDARD_REQUIRED ON
                                CXX_EXTENSIONS OFF)

target_compile_options(SU3_Lattice PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall -Wextra -pedantic>
                                         $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)