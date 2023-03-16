# Add the library target for jots
add_library(jotsLib conduction_operator.cpp
			  config_file.cpp
			  boundary_condition.cpp
			  thermal_diffusivity.cpp
			  jots_driver.cpp)

# Link MFEM as conduction_operator needs it
target_link_libraries(jotsLib PRIVATE mfem)

# Link Boost headers library as config_file.cpp needs it
target_link_libraries(jotsLib PRIVATE Boost::headers)