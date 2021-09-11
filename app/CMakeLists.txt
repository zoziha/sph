set(prog
  sph.f
)

add_executable(
  "${PROJECT_NAME}-exe"
  "${prog}"
)

target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)