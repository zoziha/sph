cmake_minimum_required(VERSION 3.14)

project(
  "SPH"
  LANGUAGES "Fortran"
  VERSION "0.0.1"
)

# Collect source of the project
add_subdirectory("src")

# Collect source of the application
add_subdirectory("app")

# Include test suite
#enable_testing()
#add_subdirectory("test")