cmake_minimum_required(VERSION 3.15)
# Set CMP0167 to NEW to use the modern FindBoost behavior
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
project(pyjeo LANGUAGES CXX)

# Must be at the start of the file
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 1. RPATH setup (Crucial for the "Build" case, ignored for the "System" case)
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# 2. Logic to find or build miallib
find_package(miallib QUIET)
if(miallib_FOUND)
    message(STATUS "Found system miallib")
else()
    message(STATUS "miallib not found, building from source")
    add_subdirectory(external/miallib)
endif()

# 3. Logic to find or build jiplib
find_package(jiplib QUIET)
if(jiplib_FOUND)
    message(STATUS "Found system jiplib")
else()
    message(STATUS "jiplib not found, building from source")
    add_subdirectory(external/jiplib)
endif()

# 1. Find the actual library file path
#find_library(JSONCPP_LIB_PATH NAMES jsoncpp PATHS ${JSONCPP_LIBRARY_DIRS})

# 4. Define the Python Extension (SWIG)
find_package(SWIG REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
include(${SWIG_USE_FILE})

# Ensure we can find headers regardless of source
target_include_directories(jiplib_python PRIVATE
    Python::Module
    #${PYTHON_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/jiplib/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/miallib/src>
    #external/jiplib/src
    #external/miallib/src
)

# 5. Smart Linking
# Link against the targets (CMake handles whether they are system or local)
target_link_libraries(jiplib_python PRIVATE jiplib miallib Python::Module)
#Ensure jiplib is compiled with a standard high enough for JsonCpp
target_compile_features(jiplib PRIVATE cxx_std_11)

# 6. Installation
# If we built them locally, we install them into the wheel root
if(NOT miallib_FOUND)
    install(TARGETS miallib DESTINATION .)
endif()
if(NOT jiplib_FOUND)
    install(TARGETS jiplib DESTINATION .)
endif()

install(TARGETS jiplib_python DESTINATION .)
install(FILES $<TARGET_PROPERTY:jiplib_python,SWIG_SUPPORT_FILES_DIRECTORY>/jiplib.py DESTINATION .)
install(DIRECTORY pyjeo/ DESTINATION pyjeo)
