add_library(rbf)

target_sources(rbf
  PUBLIC
    rbf.f90
)

target_link_libraries(rbf
  PUBLIC
    core
)