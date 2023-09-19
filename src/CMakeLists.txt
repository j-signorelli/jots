# Add the library target for jots
add_library(jotsLib conduction_operator.cpp
			  config_file.cpp
			  boundary_condition.cpp
			  material_property.cpp
			  precice_adapter.cpp
			  output_manager.cpp
			  jots_driver.cpp)

target_include_directories(jotsLib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(jotsLib PRIVATE mfem Boost::headers precice::precice)