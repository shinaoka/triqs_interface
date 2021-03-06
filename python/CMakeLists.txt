# where will the python end up in triqs?
set(python_destination pytriqs/applications/impurity_solvers/alps_cthyb)

# site_customize for build
set(package_name "pytriqs.applications.impurity_solvers")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# make the build_xxx, install python files...
triqs_prepare_local_pytriqs(${python_destination})

# Build the python module for the solver
triqs_python_extension(alps_cthyb ${python_destination})
target_link_libraries(alps_cthyb alps_cthyb_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(alps_cthyb)
install(TARGETS alps_cthyb DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
