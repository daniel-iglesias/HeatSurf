IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(VTK_FOUND)

FIND_PACKAGE(OpenGL REQUIRED)

FIND_PACKAGE(GSL REQUIRED)

# Make sure the compiler can find include files from our directories.
INCLUDE_DIRECTORIES(
  ${HeatSurf_SOURCE_DIR}/src 
  ${VTK_INCLUDE_DIR}
  ${GSL_INCLUDE_DIRS}
)

SET(HeatSurfSrcs 
  heatsurf.cpp
  cone.cpp
  cone2.cpp
  cone3.cpp
  cylinder.cpp
  cornerCylinder.cpp
  edgeCylinder.cpp
  element.cpp
  geometry.cpp
  node.cpp
  ogive.cpp
  particle.cpp
  ring.cpp
  simulation.cpp
  twoplates.cpp
)

ADD_EXECUTABLE(heatsurf
  ${HeatSurfSrcs}
)

INSTALL(TARGETS heatsurf DESTINATION bin)

link_directories(${VTK_DIR})

# Link the executable to the VTK libraries.
TARGET_LINK_LIBRARIES(heatsurf
  ${VTK_LIBRARIES}
  ${GSL_LIBRARIES}
)
