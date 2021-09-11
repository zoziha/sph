set(SRC
  art_heat.f
  art_visc.f
  av_vel.f
  density.f
  direct_find.f
  eos.f
  external_force.f
  grid_geom.f
  hsml.f
  init_grid.f
  input.f
  internal_force.f
  kernel.f
  link_list.f
  output.f
  single_step.f
  time_integration.f
  virt_part.f
  viscosity.f
  time_elapsed.f90
  time_print.f90
)

add_library(${PROJECT_NAME}-lib ${SRC})