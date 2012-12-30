set (thermal_balance_SOURCES
    therm_bal.f90
    )

add_executable (${PROJECT_NAME} ${thermal_balance_SOURCES})

find_package(LAPACK REQUIRED)
set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
