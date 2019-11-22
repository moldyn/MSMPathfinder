cmake_minimum_required (VERSION 2.8)
project (MSMPathfinder)

# do we want static libraries?
# When STATIC_LINKING is set, cmake looks for libraries called *.a.
if(STATIC_LINKING)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif(STATIC_LINKING)

# Gnu compiler setup
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  message("using Gnu compiler")
  # use C++11 standard
  set (CMAKE_CXX_FLAGS "-std=c++11")
  # profiling
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
  # optimization
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
  # parallelization
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  # warnings
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-declarations")
  # for gcc linting
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsyntax-only")


  if (NATIVE_COMPILATION)
    message("native compilation: this code will ONLY run on THIS machine")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
  set (DC_MEM_ALIGNMENT 32)
  message("using standard CPU acceleration with memory aligned at ${DC_MEM_ALIGNMENT} bytes")
else()
  message(SEND_ERROR "currently no other compiler than the GCC is supported. sorry.")
endif()



#configure_file(config.hpp.cmake.in ${CMAKE_BINARY_DIR}/generated/config.hpp)
include_directories(${CMAKE_BINARY_DIR}/generated/)

find_package (Boost 1.49 COMPONENTS program_options system filesystem REQUIRED)
include_directories (${Boost_INCLUDE_DIR})



set(CLUSTERING_SRCS src/main.cpp
                    src/paths.cpp
                    src/pathfinder.cpp
                    src/transition_matrix.cpp
                    src/progressbar.cpp)

set(CLUSTERING_LIBS ${Boost_LIBRARIES})

set (PROGNAME MSMPathfinder)
add_executable(${PROGNAME} ${CLUSTERING_SRCS})
target_link_libraries(${PROGNAME} ${CLUSTERING_LIBS})


# set -static, when STATIC_LINKING is TRUE and set LINK_SEARCH_END_STATIC
# to remove the additional -bdynamic from the linker line.
if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set_target_properties(${PROGNAME} PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(STATIC_LINKING)


install(TARGETS ${PROGNAME} RUNTIME DESTINATION .)
