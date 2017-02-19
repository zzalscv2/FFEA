
set(EIGEN_INSTALL_DIR ${EIGEN3_HOME})

include(ExternalProject)
ExternalProject_Add(getEigen
                    URL https://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2
                    URL_HASH SHA256=3e1fa6e8c45635938193f84fee6c35a87fac26ee7c39c68c230e5080c4a8fe98
                    BUILD_IN_SOURCE 0
                    CONFIGURE_COMMAND ${CMAKE_COMMAND} ../getEigen -DCMAKE_INSTALL_PREFIX=${EIGEN_INSTALL_DIR} 
                    BUILD_COMMAND "" 
                    INSTALL_COMMAND make install 
                   )

set(EIGEN3_INCLUDE_DIR ${EIGEN3_HOME}/include/eigen3)
