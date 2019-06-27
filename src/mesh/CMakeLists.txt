# mesh/CMakeLists.txt

ADD_LIBRARY(mesh
  SurfaceMesh.hpp
  utils.hpp
  Mesh.hpp
  ShapeID.hpp
  cell_iterator.hpp
  face_iterator.hpp
  mesh_methods.hpp
  surface_mesh_methods.hpp
  MetisData.hpp
  SupportRegion.hpp
  # implementation
  ShapeID.cpp
  Mesh.cpp
  cell_iterator.cpp
  face_iterator.cpp
  mesh_methods.cpp
  surface_mesh_methods.cpp
  MetisData.cpp
  SupportRegion.cpp
)

SET_TARGET_PROPERTIES (
    mesh
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC
	${CMAKE_SOURCE_DIR}/src
	${CMAKE_SOURCE_DIR}/src/uint256
  ${Boost_INCLUDE_DIR}
  ${METIS_INCLUDE_DIRS}
  ${CGAL_INCLUDE_DIRS}
)





if(0) #Boost_FOUND)
  TARGET_LINK_LIBRARIES(mesh angem ${Boost_LIBRARIES} ${METIS_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(mesh angem uint256 ${METIS_LIBRARIES} ${CGAL_LIBRARIES})
endif()
