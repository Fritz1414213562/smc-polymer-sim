target_sources(opencm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ExcludedVolumeForceFieldGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PolynomialRepulsiveForceFieldGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HarmonicAngleForceFieldGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HarmonicBondForceFieldGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MonteCarloAnisotropicBarostatGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PullingForceFieldGenerator.cpp
    )
