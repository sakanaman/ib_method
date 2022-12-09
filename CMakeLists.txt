cmake_minimum_required(VERSION 3.15)

project(fluid_experiment)

enable_language(Fortran)

find_package(OpenMP REQUIRED)

set(CMAKE_FORTRAN_FLAGS "-O3")

add_executable(single_circle
    src/interface.f90
    src/global.f90
    src/discrete_space.f90
    src/solve_ellipse.f90
    src/prepare_ellipse.f90
    src/ibm.f90
    src/initial_condition.f90
    src/boundary.f90
    src/util.f90
    src/logger.f90
    main.f90)

target_link_libraries(single_circle OpenMP::OpenMP_Fortran)
