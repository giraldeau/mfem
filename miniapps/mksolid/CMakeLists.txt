# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(mksolidlib STATIC
    spherecoeff.cpp
    tetgen.cxx
    predicates.cxx
)
target_link_libraries(mksolidlib PUBLIC mfem mfem-common Eigen3::Eigen)
target_compile_definitions(mksolidlib PUBLIC -DTETLIBRARY -DMC_IMPLEM_ENABLE)

add_mfem_miniapp(mksolid
    MAIN mksolid.cpp LIBRARIES mksolidlib)

add_mfem_miniapp(scratch
    MAIN scratch.cpp LIBRARIES mksolidlib)

add_mfem_miniapp(impsphere
    MAIN impsphere.cpp
    ${MFEM_MINIAPPS_COMMON_HEADERS}
    LIBRARIES mksolidlib)

add_mfem_miniapp(mcsphere
    MAIN mcsphere.cpp
    ${MFEM_MINIAPPS_COMMON_HEADERS}
    LIBRARIES mksolidlib)

add_mfem_miniapp(triwedge
    MAIN triwedge.cpp
    ${MFEM_MINIAPPS_COMMON_HEADERS}
    LIBRARIES mksolidlib)

add_mfem_miniapp(surfmesh
    MAIN surfmesh.cpp
    ${MFEM_MINIAPPS_COMMON_HEADERS}
    LIBRARIES mksolidlib)
