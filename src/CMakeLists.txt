
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common sdsl malloc_count)

add_executable(rlbwt rlbwt.cpp)
target_link_libraries(rlbwt common sdsl malloc_count)

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

add_executable(query_mum_phinder mum-phinder.cpp)
 target_link_libraries(query_mum_phinder common sdsl divsufsort divsufsort64 malloc_count ri)
 target_include_directories(query_mum_phinder PUBLIC    "../include/ms" 
                                         "../include/common" 
                                         "${shaped_slp_SOURCE_DIR}" 
                                         "${FOLCA_SOURCE_DIR}" 
                                         "${SUX_SOURCE_DIR}/function" 
                                         "${SUX_SOURCE_DIR}/support"
                                         )
 target_compile_options(query_mum_phinder PUBLIC "-std=c++17")

add_executable(build_mum_phinder build_mum-phinder.cpp)
target_link_libraries(build_mum_phinder common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(build_mum_phinder PUBLIC    "../include/ms" 
                                        "../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(build_mum_phinder PUBLIC "-std=c++17")