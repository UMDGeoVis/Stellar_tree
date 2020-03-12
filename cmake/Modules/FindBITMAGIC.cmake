set(BM_PREFIX "" CACHE PATH "path ")


find_path(BM_INCLUDE_DIR bm/bm.h
    PATHS ${BM_PREFIX}/include /usr/include /usr/local/include ${BITMAGIC_DIR}/include)

#find_library(GMP_LIBRARY NAMES gmp libgmp 
#    PATHS ${BM_PREFIX}/lib /usr/lib /usr/local/lib)

if(BM_INCLUDE_DIR)
    #get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
    set(BM_FOUND TRUE)
endif()

if(BM_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found BitMagic: ${BM_INCLUDE_DIR}")
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find BitMagic")
   endif()
endif()
