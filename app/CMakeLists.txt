set(prog
  sph.f90
)

add_executable(
  "${PROJECT_NAME}-exe"
  "${prog}"
)

set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)

target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)