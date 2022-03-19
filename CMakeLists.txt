cmake_minimum_required(VERSION 3.16)

project(Bode VERSION 0.0.1)

find_package(ROOT)
find_library(ERR_A_LIB ErrorAnalysis)
find_path(ERR_A_PATH ErrorAnalysis.h PATHS /usr/local)

find_library(LAB_LIB LabTools)
find_path(LAB_PATH
    LabTools.h
    LabPlot.h
    LabStat.h
    ATLASStyle.h
    Logger.h 
    PATHS /usr/local)

set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})

set(BODEINC
    BodeDataSim/SimEngine.h
    Bode/Analysis.h)
    
set(BODESRC
    src/Analysis.cpp)

add_compile_options(-I${ROOT_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(Bode SHARED ${BODESRC})

target_link_libraries(Bode ${ROOT_LIBRARIES} ${ERR_A_LIB} ${LAB_LIB})
target_include_directories(Bode PUBLIC ${ERR_A_PATH} ${LAB_PATH})

install(TARGETS Bode DESTINATION /usr/local/lib/)
install(FILES ${BODEINC} DESTINATION include/Bode)


