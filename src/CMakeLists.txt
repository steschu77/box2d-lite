project(box2d-lite LANGUAGES CXX)

set(BOX2D_SOURCE_FILES
	Arbiter.cpp
    Arbiter.h
	Body.cpp
    Body.h
	Collide.cpp
	Joint.cpp
    Joint.h
    Matrix2x2.cpp
    Vector2.cpp
    VecMath.h
	World.cpp
    World.h
)

add_library(box2d-lite STATIC ${BOX2D_SOURCE_FILES})
target_include_directories(box2d-lite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
