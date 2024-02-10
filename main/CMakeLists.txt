add_executable(interpolator main.f90)
target_link_libraries(interpolator
  PRIVATE
    core
)