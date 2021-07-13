# Generated on 2020-2-20 16:14:09 by version 0.9.0 of mulle-xcode-to-cmake
# Command line:
#    mulle-xcode-to-cmake -b -l c++ export ../IntervalAssignment-xcode/IntervalAssignment-xcode.xcodeproj


cmake_minimum_required (VERSION 3.4)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project( test-IntervalAssignment CXX)


##
## test-IntervalAssignment Files
##

# uncomment this for mulle-objc to search libraries first
# set( CMAKE_FIND_FRAMEWORK "LAST")

set( SOURCES
../Source/CpuTimer.cpp
../Source/IA.cpp
../Source/IAQueue.cpp
../Source/IAResultImplementation.cpp
../Source/IncrementalIntervalAssignment.cpp
../Source/test.cpp
)

set( STATIC_DEPENDENCIES
)

##
## test-IntervalAssignment
##

add_executable( test-IntervalAssignment
${SOURCES}
${PUBLIC_HEADERS}
${PROJECT_HEADERS}
${PRIVATE_HEADERS}
${RESOURCES}
)

target_link_libraries( test-IntervalAssignment
${STATIC_DEPENDENCIES}
${DEPENDENCIES}
)
