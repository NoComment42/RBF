cmake_minimum_required(VERSION 3.26)

project(interplation
  LANGUAGES
    Fortran
)

# Set cmake paths up in my preferred way. This makes it easier to find/link against
# multiple modules/libraries.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

set(ENABLE_LAPACK95 True)
find_package(MKL CONFIG REQUIRED)

add_subdirectory(core)
add_subdirectory(rbf)
add_subdirectory(main)