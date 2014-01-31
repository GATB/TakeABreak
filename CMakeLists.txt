project(TakeABreak)

cmake_minimum_required (VERSION 2.6)

################################################################################
# GATB-CORE
################################################################################

# we depend on gatb-core; here, we define where to find all the required material
add_subdirectory (../../thirdparty/gatb-core/gatb-core  "${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core")


################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we give the headers directories from : 
#       - from project source 
#       - from GATB-CORE source
#       - from dsk source 
include_directories (include ${gatb-core-includes})

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

#file (GLOB_RECURSE  TestDBGWalkerFiles  src/* tests/test_DBGWalker_texas_ranger.cpp)

#file (GLOB_RECURSE  TestClavetteFiles  src/* tests/testClavette.cpp)

# we define the artifact to be built: the project binary
add_executable        (${CMAKE_PROJECT_NAME}  src/${CMAKE_PROJECT_NAME}.cpp ${ProjectFiles})
#add_executable        (testDBGWalker ${TestDBGWalkerFiles})
#add_executable        (testClavette ${TestClavetteFiles})

# we define which libraries to be linked with project binary
target_link_libraries (${CMAKE_PROJECT_NAME}  ${gatb-core-libraries})
#target_link_libraries (testDBGWalker  ${gatb-core-libraries})
#target_link_libraries (testClavette  ${gatb-core-libraries})

