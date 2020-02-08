cmake_minimum_required(VERSION 3.8)

project(box2d-lite LANGUAGES CXX)

set(BOX2D_SOURCE_FILES
	src/Arbiter.cpp
    src/Arbiter.h
	src/Body.cpp
    src/Body.h
	src/Collide.cpp
	src/Joint.cpp
    src/Joint.h
    src/Matrix2x2.cpp
    src/Rotate2.cpp
    src/Vector2.cpp
    src/VecMath.h
	src/World.cpp
    src/World.h
)

add_library(box2d-lite STATIC ${BOX2D_SOURCE_FILES})
target_include_directories(box2d-lite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

option(BOX2D_BUILD_SAMPLES "Build the box2d-lite sample program" ON)

if (BOX2D_BUILD_SAMPLES)

	add_subdirectory(extern/glad)
	add_subdirectory(extern/glfw)
	add_subdirectory(extern/imgui)
	add_subdirectory(samples)

	# default startup project for Visual Studio
	if (MSVC)
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT samples)
	endif()

endif()
