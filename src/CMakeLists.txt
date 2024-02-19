# Add the library target for jots
add_library(jotsLib 
                jots_iterator.cpp
                linear_conduction_operator.cpp
                nl_conduction_operator.cpp
                steady_conduction_operator.cpp
                config_file.cpp
                boundary_condition.cpp
                material_property.cpp
                jots_precice.cpp
                output_manager.cpp
                jots_driver.cpp
                jots_nlfis.cpp)

target_include_directories(jotsLib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


# Link dependencies and any additional compiler flags
set(dependencies mfem Boost::headers precice::precice)
target_link_libraries(jotsLib PUBLIC ${dependencies} compiler_flag_requirements)

# Install
set(installable_libs jotsLib compiler_flag_requirements)
install(TARGETS ${installable_libs} DESTINATION lib)
install(FILES "jots_driver.hpp" DESTINATION include)