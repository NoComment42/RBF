add_library(rbf)

target_sources(rbf
  PUBLIC
    rbf.f90
)

target_link_libraries(rbf
  PRIVATE
    core
)

if( ${CMAKE_Fortran_COMPILER_ID} STREQUAL IntelLLVM )
  message("Using IntelLLVM compiler")
  target_include_directories(rbf PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_compile_options(rbf PUBLIC -fpp $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_link_libraries(rbf
    PUBLIC
      $<LINK_ONLY:MKL::MKL>
  )
endif()