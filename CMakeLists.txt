#----------------------------------------------------------------------------
# Setup the project
#

## Set the default build type
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Set build type (RelWithDebInfo is default)")

set(DEV_WARNING "cmake build of SBS-Offline is still in development. Please test and report issues")

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(sbs)

# Setup GNUInstallDirs for installation:
#----------------------------------------------------------------------------
include(GNUInstallDirs)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------

## Disable in-source builds!
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "cmake cannot build SBS-offline in the root directory: "
    "${PROJECT_SOURCE_DIR}\nsince it conflicts with our old Makefile."
    " Create a build directory anywhere, change to it and re-run cmake."
    )
endif()

##----------------------------------------------------------------------------
## Set the sources which have a corresponding .h file here
set(sources MPDModule.cxx MPDModuleVMEv4.cxx SBSBigBite.cxx
  SBSGEMTrackerBase.cxx SBSGEMSpectrometerTracker.cxx SBSGEMModule.cxx
  SBSBBShower.cxx SBSBBTotalShower.cxx 
  SBSCDet.cxx SBSScintHit.cxx SBSScintPMT.cxx 
  SBSTimingHodoscope.cxx SBSScintBar.cxx SBSTdcHit.cxx SBSAdcHit.cxx
  SBSScintPartialHit.cxx SBSCherenkovDetector.cxx SBSCherenkov_ClusterList.cxx
  SBSGRINCH.cxx
  SBSScintPlane.cxx SBSEArm.cxx
  SBSHCal.cxx SBSDecodeF1TDCModule.cxx SBSCalorimeter.cxx SBSGenericDetector.cxx
  SBSData.cxx SBSElement.cxx
  SBSCalorimeterCluster.cxx SBSSimDataDecoder.cxx 
  SBSSimDecoder.cxx SBSSimADC.cxx SBSSimTDC.cxx
  SBSHCalLEDModule.cxx SBSManager.cxx
  SBSSimFile.cxx SBSSimEvent.cxx
  SBSRPBeamSideHodo.cxx SBSRPFarSideHodo.cxx SBSCHAnalyzer.cxx
  SBSTimingHodoscopePMT.cxx SBSTimingHodoscopeBar.cxx SBSTimingHodoscopeCluster.cxx
  SBSBPM.cxx SBSRaster.cxx SBSRasteredBeam.cxx 
  LHRSScalerEvtHandler.cxx SBSScalerEvtHandler.cxx
  SBSScalerHelicity.cxx SBSScalerHelicityReader.cxx
  gmn_tree_digitized.cxx genrp_tree_digitized.cxx VETROC.cxx
  )
#gmn_tree_digitized.C
#g4sbs_data.cxx g4sbs_tree.cxx
#gmn_dig_tree.cxx  SBSSimMPD.cxx SBSBBShowerCluster.cxx SBSShowerBlock.cxx
#SBSGRINCH.cxx SBSGRINCH_ClusterList.cxx

## Make corresponding header files
string(REGEX REPLACE .cxx .h headers "${sources}")
## Append a few extra header files in the repository.  (But do not append the
## sbs_LinkDef.h header here, that is further down this cmake config.
list(APPEND headers "DebugDef.h" "Helper.h")

## Setup default variables (copied from Makefile)
#FIXME: should be configure'd:
set(EXTRAWARN ON CACHE BOOL "g++ Wextra warn options if supported")

##----------------------------------------------------------------------------
set(VERBOSE ON CACHE BOOL "Compile extra code for printing verbose messages")
set(TESTCODE ON CACHE BOOL "Compile extra diagnostic code (extra computations and global variables")
set(MCDATA ON CACHE BOOL "Compile support code for MC input data")
set(CXXMAXERRORS 0 CACHE STRING "Maximum number of allowed errors before cxx stops")
list(APPEND SBSEXTRADEF_LIST VERBOSE TESTCODE MCDATA)

#----------------------------------------------------------------------------
# Find ROOT 
#ROOT path
# commented out the following line so that this procedure will evolve with ROOT:
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
##----------------------------------------------------------------------------
## Setup the PODD Analyzer
find_package(Podd 1.7 REQUIRED)

##----------------------------------------------------------------------------
## Check for SStream support
include(TestForSSTREAM)
if(CMAKE_HAS_ANSI_STRING_STREAM)
  set(HAS_SSTREAM ON)
  list(APPEND SBSEXTRADEF_LIST HAS_SSTREAM)
endif()

##----------------------------------------------------------------------------
## OS Type (not sure what exactly uses this, but since it's in the Makefile
## I'll just copy it here too...
set(SBS_LINKFLAGS "")
if(UNIX) ## Can be anything from Linux to Mac OS X
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")
  if(APPLE)
    set(MACVERS ON)
    list(APPEND SBSEXTRADEF_LIST MACVERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
    if(EXTRAWARN)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
    endif()
    set(SBS_LINKFLAGS "${SBS_LINKFLAGS} -Wl,-undefined,dynamic_lookup")
  else()
    set(LINUXVERS ON)
    list(APPEND SBSEXTRADEF_LIST LINUXVERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
    if(EXTRAWARN AND CMAKE_COMPILER_IS_GNUCXX AND NOT
        CMAKE_CXX_COMPILER_VERSION LESS 4)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
    endif()
    if(CXXMAXERRORS GREATER 0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=${CXXMAXERRORS}")
    endif()
  endif()
endif() ## I guess no Windows support?



##----------------------------------------------------------------------------
## Copy over some flags we have set in the Makefile
if(VERBOSE)
  list(APPEND VERBOSE)
endif()



include_directories(${ROOT_INCLUDE_DIR} ${Podd_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

## Generate the ROOT dictionary
set(PROJECT_DICTIONARY ${PROJECT_NAME}Dict)
ROOT_GENERATE_DICTIONARY( ${PROJECT_DICTIONARY} ${headers} OPTIONS -I${PROJECT_SOURCE_DIR} MODULE ${PROJECT_NAME} LINKDEF ${PROJECT_NAME}_LinkDef.h )

## Create the library
add_library(${PROJECT_NAME} SHARED ${sources} ${headers} ${PROJECT_DICTIONARY}.cxx )
set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY LINK_FLAGS " ${SBS_LINKFLAGS}")
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})
#message(FATAL_ERROR "${PROJECT_NAME} PRIVATE ${SBSEXTRADEFINES}")
## Get a list of preprocessor flags to pass to the targets
foreach(_def ${SBSEXTRADEF_LIST})
  if(${${_def}})
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${_def})
  endif()
endforeach()

## default SBS_REPLAY_PATH to 

# cmake_path wasn't introduced until cmake 3.20, not yet available on ifarm by default
#cmake_path(SET DEFAULT_SBS_REPLAY_PATH "${PROJECT_SOURCE_DIR}/../SBS-replay")

option( GEN "Use GEN-replay instead of SBS-replay as default replay directory" false )

set(REPLAYDIR SBS-replay)
if( GEN )
    set(REPLAYDIR GEN-replay)
endif()

if( DEFINED ENV{SBS_REPLAY} )
    set( DEFAULT_SBS_REPLAY_PATH $ENV{SBS_REPLAY} )
else() # look in directory 
   get_filename_component( DEFAULT_SBS_REPLAY_PATH ../${REPLAYDIR} ABSOLUTE BASE_DIR ${PROJECT_SOURCE_DIR})
endif()		   

set(SBS_REPLAY_PATH ${DEFAULT_SBS_REPLAY_PATH} CACHE PATH "Path to ${REPLAYDIR} installation")

# check if sbs_replay_path set by default from environment matches SBS-replay or GEN-replay, adjust if needed:
if( NOT ${SBS_REPLAY_PATH} MATCHES "/${REPLAYDIR}$" )
  if( ${SBS_REPLAY_PATH} MATCHES "/SBS-replay$" )
    message( WARNING "default replay directory ${SBS_REPLAY_PATH} from environment matches SBS-replay but GEN-replay option was given on the command-line, switching to SBS-replay" )
    set(REPLAYDIR SBS-replay)
  elseif( ${SBS_REPLAY_PATH} MATCHES "/GEN-replay$" )
    message( WARNING "default replay directory ${SBS_REPLAY_PATH} from environment matches GEN-replay but SBS-replay option was given on the command-line, switching to GEN-replay" )
    set(REPLAYDIR GEN-replay)
  else() 
    message( WARNING "default replay directory ${SBS_REPLAY_PATH} from environment doesn't match a supported configuration, appending ${REPLAYDIR}, YMMV..." )
  endif()
# 
  string(APPEND SBS_REPLAY_PATH "/${REPLAYDIR}")
endif()

if(NOT EXISTS ${SBS_REPLAY_PATH})
  if( EXISTS ${DEFAULT_SBS_REPLAY_PATH} )
    message( WARNING "${REPLAYDIR} not found in ${SBS_REPLAY_PATH}, reverting to replay installation found in ${DEFAULT_SBS_REPLAY_PATH}" )
    set(SBS_REPLAY_PATH ${DEFAULT_SBS_REPLAY_PATH} CACHE PATH "Path to ${REPLAYDIR} installation" FORCE)
  elseif(${DEFAULT_SBS_REPLAY_PATH} STREQUAL ${SBS_REPLAY_PATH})
    message(WARNING "${REPLAYDIR} not found in default location ${DEFAULT_SBS_REPLAY_PATH}. You can still build SBS-offline, but the \"replay\" and \"scripts\" directories under ${REPLAYDIR} won't be automatically added to ROOT's macro path in the .rootrc file. To fix this, set environment variable SBS_REPLAY to point to your ${REPLAYDIR} installation and/or re-run cmake with -DSBS_REPLAY_PATH=/path/to/your/custom/${REPLAYDIR}/")
  else()
    message(WARNING "${REPLAYDIR} not found in either its default location ${DEFAULT_SBS_REPLAY_PATH} or user-specified location ${SBS_REPLAY_PATH}. You can still build SBS-offline, but the \"replay\" and \"scripts\" directories under ${REPLAYDIR} won't be automatically added to ROOT's macro path in the .rootrc file. To fix this, set environment variable SBS_REPLAY to point to your ${REPLAYDIR} installation and/or re-run cmake with -DSBS_REPLAY_PATH=/path/to/your/custom/${REPLAYDIR}/")
  endif()	   
endif()	

if( EXISTS ${SBS_REPLAY_PATH} )
    set(SBS_MACRO_PATH ".:${SBS_REPLAY_PATH}/replay:${SBS_REPLAY_PATH}/scripts:${SBS_REPLAY_PATH}/onlineGUIconfig: ${SBS_REPLAY_PATH}/onlineGUIconfig/scripts")
else()
    set(SBS_MACRO_PATH ".")
endif()

configure_file( ${PROJECT_SOURCE_DIR}/sbsenv.csh ${CMAKE_CURRENT_BINARY_DIR}/sbsenv.csh @ONLY)
configure_file( ${PROJECT_SOURCE_DIR}/sbsenv.sh ${CMAKE_CURRENT_BINARY_DIR}/sbsenv.sh)

configure_file( ${PROJECT_SOURCE_DIR}/SBSoffline.rootrc ${CMAKE_CURRENT_BINARY_DIR}/.rootrc )
configure_file( ${PROJECT_SOURCE_DIR}/rootlogon_install.C ${PROJECT_BINARY_DIR}/rootlogon.C )

## Time to disable this warning. The cmake build works...
#message(AUTHOR_WARNING ${DEV_WARNING})

#install the library under ${CMAKE_INSTALL_PREFIX}/lib
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
#install the relevant directories: DB, replay, scripts in the top-level installation folder (these have been moved to ${REPLAYDIR} repository)
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/DB DESTINATION ${CMAKE_INSTALL_PREFIX})
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/DB_UVA_EEL_COSMIC_JAN2021 DESTINATION ${CMAKE_INSTALL_PREFIX})
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/replay DESTINATION ${CMAKE_INSTALL_PREFIX})
#INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX})
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbsenv.sh DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbsenv.csh DESTINATION ${CMAKE_INSTALL_BINDIR})
#install all headers under cmake_install_prefix/include
INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install the ROOT dictionary files: eventually we should learn how to do this in a smarter way. Put the .cxx file under "include"
# and put the _rdict.pcm and .rootmap files under "lib". This seems to work. Don't ask me why
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_DICTIONARY}.cxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/rootlogon.C DESTINATION ${CMAKE_INSTALL_SYSCONFDIR} )

#this creates run_replay_here folder under install prefix
INSTALL(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/run_replay_here)

# put the .rootrc file under run_replay_here
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/.rootrc DESTINATION ${CMAKE_INSTALL_PREFIX}/run_replay_here)

#option(SET_MCFLAG "Set the global MCDATA flag" OFF)
#if(SET_MCFLAG)
#endif()
