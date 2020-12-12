# ---------------------------------------------------------------------------------------
# IDE support for headers
# ---------------------------------------------------------------------------------------

# spd *h files
set(SPD_HEADERS_DIR "${CMAKE_CURRENT_LIST_DIR}/../spdlog/include")
file(GLOB_RECURSE SPD_HEADERS "${SPD_HEADERS_DIR}/*.h")
set(SPD_ALL_HEADERS "${SPD_HEADERS}")

set(HEADER_DIR ${SPD_HEADERS_DIR})
set(HDADER_NAME "spdlog")
set(HEADERS ${SPD_HEADERS})
foreach(_source IN ITEMS ${HEADERS})
    get_filename_component(_source_path "${_source}" PATH)
    file(RELATIVE_PATH _source_path_rel "${HEADER_DIR}" "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
    source_group("Header Files\\${HDADER_NAME}\\${_group_path}" FILES "${_source}")
endforeach()


# tbb *h files
set(TBB_HEADERS_DIR "${CMAKE_CURRENT_LIST_DIR}/../tbb/include")
file(GLOB_RECURSE TBB_HEADERS "${TBB_HEADERS_DIR}/*.h")
set(TBB_ALL_HEADERS "${TBB_HEADERS}")

set(HEADER_DIR ${TBB_HEADERS_DIR})
set(HDADER_NAME "tbb")
set(HEADERS ${TBB_HEADERS})
foreach(_source IN ITEMS ${HEADERS})
    get_filename_component(_source_path "${_source}" PATH)
    file(RELATIVE_PATH _source_path_rel "${HEADER_DIR}" "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
    source_group("Header Files\\${HDADER_NAME}\\${_group_path}" FILES "${_source}")
endforeach()


# tinyobjloader *h files
set(TINY_OBJ_LOADER_HEADERS_DIR "${CMAKE_CURRENT_LIST_DIR}/../tinyobjloader")
file(GLOB_RECURSE TINY_OBJ_LOADER_HEADERS "${TINY_OBJ_LOADER_HEADERS_DIR}/*.h")
set(TINY_OBJ_LOADER_ALL_HEADERS "${TINY_OBJ_LOADER_HEADERS}")

set(HEADER_DIR ${TINY_OBJ_LOADER_HEADERS_DIR})
set(HDADER_NAME "tinyobjloader")
set(HEADERS ${TINY_OBJ_LOADER_HEADERS})
foreach(_source IN ITEMS ${HEADERS})
    get_filename_component(_source_path "${_source}" PATH)
    file(RELATIVE_PATH _source_path_rel "${HEADER_DIR}" "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
    source_group("Header Files\\${HDADER_NAME}\\${_group_path}" FILES "${_source}")
endforeach()


# stb *h files
set(STB_HEADERS_DIR "${CMAKE_CURRENT_LIST_DIR}/../stb")
file(GLOB_RECURSE STB_HEADERS "${STB_HEADERS_DIR}/*.h")
set(STB_ALL_HEADERS "${STB_HEADERS}")

set(HEADER_DIR ${STB_HEADERS_DIR})
set(HDADER_NAME "stb")
set(HEADERS ${STB_HEADERS})
foreach(_source IN ITEMS ${HEADERS})
    get_filename_component(_source_path "${_source}" PATH)
    file(RELATIVE_PATH _source_path_rel "${HEADER_DIR}" "${_source_path}")
    string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
    source_group("Header Files\\${HDADER_NAME}\\${_group_path}" FILES "${_source}")
endforeach()
