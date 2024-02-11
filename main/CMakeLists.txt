add_executable(interpolator)
target_sources(interpolator
  PRIVATE
    main.f90
    helper.f90
)

target_link_libraries(interpolator
  PRIVATE
    core
    rbf
)