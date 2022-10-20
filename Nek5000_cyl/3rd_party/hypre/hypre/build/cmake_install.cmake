# Install script for directory: /home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/build/libHYPRE.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/build/HYPRE_config.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/HYPREf.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/HYPRE.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/HYPRE_utilities.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/_hypre_utilities.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/hypre_hopscotch_hash.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/HYPRE_error_f.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/fortran.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/utilities/fortran_matrix.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/csr_matmultivec.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/interpreter.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/multivector.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/par_csr_matmultivec.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/par_csr_pmvcomm.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/par_multivector.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/seq_multivector.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/multivector/temp_multivector.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/krylov/HYPRE_krylov.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/krylov/HYPRE_lobpcg.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/krylov/HYPRE_MatvecFunctions.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/krylov/krylov.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/krylov/lobpcg.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/seq_mv/HYPRE_seq_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/seq_mv/seq_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_mv/HYPRE_parcsr_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_mv/_hypre_parcsr_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_block_mv/par_csr_block_matrix.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_block_mv/csr_block_matrix.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/distributed_matrix/distributed_matrix.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/IJ_mv/HYPRE_IJ_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/IJ_mv/_hypre_IJ_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/matrix_matrix/HYPRE_matrix_matrix_protos.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/distributed_ls/pilut/HYPRE_DistributedMatrixPilutSolver_protos.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/distributed_ls/pilut/HYPRE_DistributedMatrixPilutSolver_types.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_ls/HYPRE_parcsr_ls.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/parcsr_ls/_hypre_parcsr_ls.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/struct_mv/HYPRE_struct_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/struct_mv/_hypre_struct_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/struct_ls/HYPRE_struct_ls.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/struct_ls/_hypre_struct_ls.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/sstruct_mv/HYPRE_sstruct_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/sstruct_mv/_hypre_sstruct_mv.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/sstruct_ls/HYPRE_sstruct_ls.h"
    "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/src/sstruct_ls/_hypre_sstruct_ls.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/prabal/workstation/git/nek_cylindrical/Nek5000/3rd_party/hypre/hypre/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
