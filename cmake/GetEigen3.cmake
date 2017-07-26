# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#


set(EIGEN_INSTALL_DIR ${EIGEN3_HOME})

include(ExternalProject)
# Just get Eigen, and forget configuring and installing the package. 
#   On some cases it was looking for Boost, causing trouble if finding an old version.   
		    # CONFIGURE_COMMAND ${CMAKE_COMMAND} ../getEigen -DCMAKE_INSTALL_PREFIX=${EIGEN_INSTALL_DIR} 
		    # INSTALL_COMMAND make install 
ExternalProject_Add(getEigen
                    URL https://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2
                    URL_HASH SHA256=dd254beb0bafc695d0f62ae1a222ff85b52dbaa3a16f76e781dce22d0d20a4a6
                    BUILD_IN_SOURCE 0
		    CONFIGURE_COMMAND ""
                    BUILD_COMMAND "" 
		    INSTALL_COMMAND ""
                   )

#set(EIGEN3_INCLUDE_DIR ${EIGEN3_HOME}/include/eigen3)
set(EIGEN3_INCLUDE_DIR ${CMAKE_BINARY_DIR}/getEigen-prefix/src/getEigen)
