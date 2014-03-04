################################################################################
 #   TakeABreak: Breakpoint detection from raw unassembled NGS reads
 #   A tool from the GATB (Genome Assembly Tool Box)
 #   Copyright (C) 2014  INRIA
 #   Authors: C.Lemaitre, P.Peterlongo, E.Drezen
 #
 #  This program is free software: you can redistribute it and/or modify
 #  it under the terms of the GNU Affero General Public License as
 #  published by the Free Software Foundation, either version 3 of the
 #  License, or (at your option) any later version.
 #
 #  This program is distributed in the hope that it will be useful,
 #  but WITHOUT ANY WARRANTY; without even the implied warranty of
 #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #  GNU Affero General Public License for more details.
 #
 #  You should have received a copy of the GNU Affero General Public License
 #  along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

project(TakeABreak)

cmake_minimum_required (VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
LIST (APPEND CMAKE_PATH_LOOKUP "cmake;../../cmake")

FOREACH (path ${CMAKE_PATH_LOOKUP})
    IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
    ENDIF()
ENDFOREACH(path)

################################################################################
# GATB-CORE CONFIGURATION
################################################################################

# We include all the required material for gatb-core
include (GatbCore)

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
#file (GLOB_RECURSE  TestClavetteFiles   src/* tests/testClavette.cpp)

# we define the artifact to be built: the project binary
add_executable        (${PROJECT_NAME}  src/${PROJECT_NAME}.cpp ${ProjectFiles})
#add_executable        (testDBGWalker ${TestDBGWalkerFiles})
#add_executable        (testClavette  ${TestClavetteFiles})

# we define which libraries to be linked with project binary
target_link_libraries (${PROJECT_NAME}  ${gatb-core-libraries})
#target_link_libraries (testDBGWalker   ${gatb-core-libraries})
#target_link_libraries (testClavette    ${gatb-core-libraries})

################################################################################
#  PACKAGING 
################################################################################

# We chose the kind of archive we want to generate
SET (CPACK_GENERATOR            "TGZ")
SET (CPACK_SOURCE_GENERATOR     "TGZ")

# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)

# We copy other stuff 
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data          DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc           DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README        DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE       DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/TakeABreak.sh DESTINATION .)

# For the source archive, we exclude some unwanted directories.
SET (CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES} "/build/"  "/distrib/"   "/tests/"  ".gitignore$")

# To be done at the end. 
INCLUDE (CPack)

