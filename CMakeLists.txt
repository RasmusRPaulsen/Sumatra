cmake_minimum_required(VERSION 3.17)

project(Sumatra)

#======================= INCLUSION OF Qt =======================#
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_PREFIX_PATH $ENV{QTDIR})
find_package(OpenGL)
find_package(Qt6Core REQUIRED)
find_package(Qt6Gui REQUIRED)
find_package(Qt6OpenGL REQUIRED)
find_package(Qt6Xml REQUIRED)
find_package(Qt6Widgets REQUIRED)

#======================= INCLUSION OF VTK ======================#
set(VTK_DIR $ENV{VTK_DIR})
find_package(VTK
    REQUIRED
    COMPONENTS
        GUISupportQt
        IOLegacy
		IOImport
		IOImage
		IOGeometry
		IOPLY
		IOXML
		IOExport
		ImagingGeneral
		ImagingHybrid
		FiltersHybrid
		FiltersModeling
		CommonColor
		CommonCore
		CommonSystem
		CommonDataModel
		RenderingAnnotation
		InteractionStyle
		RenderingContextOpenGL2
		RenderingCore
		RenderingFreeType
		RenderingGL2PSOpenGL2
		RenderingOpenGL2)
				

#======================= INCLUSION OF GEL ======================#
# Could probably be made much better
# set(GEL_DIR $ENV{GEL_SRC_DIR})				
# FIND_PATH(GEL_SRC_DIR GEL/HMesh/Manifold.h ${GEL_DIR}/src )
# include_directories(${GEL_SRC_DIR})	

# set(GEL_BDIR $ENV{GEL_BIN_DIR})				
# FIND_LIBRARY(GEL_LIBRARY_RELEASE GEL.lib ${GEL_BDIR}/Release )
# FIND_LIBRARY(GEL_LIBRARY_DEBUG GEL.lib ${GEL_BDIR}/Debug )
				
#======================= INCLUSION OF Our Code ======================#
set(FORMS_DIR "${CMAKE_SOURCE_DIR}/forms")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(RESOURCE_DIR "${CMAKE_SOURCE_DIR}/resources")

include_directories(${FORMS_DIR})
include_directories(${INCLUDE_DIR})
include_directories(${SOURCE_DIR})
include_directories(${RESOURCE_DIR})

file(GLOB_RECURSE SOURCES
    "${FORMS_DIR}/*.ui"
    "${INCLUDE_DIR}/*.h"
    "${SOURCE_DIR}/*.cpp"
    "${SOURCE_DIR}/*.cxx"
    "${RESOURCE_DIR}/*.qrc"
)

# To enable windows icon
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    enable_language("RC")
    set (WIN32_RESOURCES ${CMAKE_SOURCE_DIR}/windows/Sumatra.rc)
endif()


# Add the forms directory to the AUTOUIC search paths
set(CMAKE_AUTOUIC_SEARCH_PATHS ${CMAKE_AUTOUIC_SEARCH_PATHS} ${FORMS_DIR})

# Add the executable
if (WIN32) 
	# Prevent a "command line is too long" failure in Windows.
	set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
    add_executable(Sumatra WIN32 ${SOURCES} ${WIN32_RESOURCES})
elseif(UNIX)
    add_executable(Sumatra ${SOURCES})
endif()

# Add the target includes
target_include_directories(Sumatra PRIVATE ${FORMS_DIR})
target_include_directories(Sumatra PRIVATE ${INCLUDE_DIR})
target_include_directories(Sumatra PRIVATE ${SOURCE_DIR})
target_include_directories(Sumatra PRIVATE ${RESOURCE_DIR})

#===================== LINKING LIBRARIES =======================#
target_link_libraries( Sumatra Qt6::OpenGL)
target_link_libraries( Sumatra Qt6::Xml)
target_link_libraries( Sumatra Qt6::Widgets)
target_link_libraries( Sumatra ${QT_LIBRARIES} ${OPENGL_LIBRARIES})
#target_link_libraries( Sumatra VTK::GUISupportQt VTK::IOLegacy VTK:IOGeometry)
target_link_libraries( Sumatra ${VTK_LIBRARIES})
# target_link_libraries(Sumatra debug ${GEL_LIBRARY_DEBUG})
# target_link_libraries(Sumatra optimized ${GEL_LIBRARY_RELEASE})

vtk_module_autoinit(
    TARGETS Sumatra
    MODULES ${VTK_LIBRARIES})
