set(SRC
  sph_kind.F90
  art_heat.f90
  art_visc.f90
  av_vel.f90
  density.f90
  direct_find.f90
  eos.f90
  external_force.f90
  grid_geom.f90
  hsml.f90
  init_grid.f90
  input.f90
  internal_force.f90
  kernel.f90
  link_list.f90
  output.f90
  single_step.f90
  time_integration.f90
  virt_part.f90
  viscosity.f90
  time_elapsed.f90
  time_print.f90
)

set(
  MODULE_OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/include"
)

add_library(${PROJECT_NAME}-lib ${SRC})

set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${MODULE_OUTPUT_DIR}"
)

target_include_directories(
  "${PROJECT_NAME}-lib"
  INTERFACE
  $<BUILD_INTERFACE:${MODULE_OUTPUT_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)