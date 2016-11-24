# The solver
add_library(alps_cthyb_c solver_core.cpp)
#target_link_libraries(alps_cthyb_c ${TRIQS_LIBRARY_ALL} ${ALPSCore_LIBRARIES} ${ALPSCoreCTHYB_LIBRARIES})
target_link_libraries(alps_cthyb_c ${TRIQS_LIBRARY_ALL} ${ALPSCore_LIBRARIES} alpscore_cthyb)
#message(${ALPSCoreCTHYB_LIBRARIES})
#message(${ALPSCore_LIBRARIES})
#message(${ALPSCoreCTHYB_INCLUDE_DIR})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
#target_link_libraries(alps_cthyb ${ALPSCore_LIBRARIES} ${ALPSCoreCTHYB_LIBRARIES})
triqs_set_rpath_for_target(alps_cthyb_c)

#configure_file(config.hpp.in config.hpp)
install(TARGETS alps_cthyb_c DESTINATION lib)

set(HEADERS_TO_INSTALL
    solve_parameters.hpp
    solver_core.hpp
    config.hpp
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.hpp DESTINATION include/alps_cthyb)
install(FILES ${HEADERS_TO_INSTALL} DESTINATION include/alps_cthyb)
