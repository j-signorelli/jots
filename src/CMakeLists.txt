# Add the library target for jots
add_library(jotsLib conduction_operator.cpp
			  config_file.cpp)

# Link MFEM as conduction_operator needs it
target_link_libraries(jotsLib PRIVATE mfem)