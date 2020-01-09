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

set (SPHEREDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/connection_identify_face/sphere_structure")
file (COPY round_lad.1.vdw DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.top DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.surf DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.stokes DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.pin DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.node DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.mat DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.lj DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.smesh DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.face DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.ele DESTINATION ${SPHEREDIR})
file (COPY round_lad.1.edge DESTINATION ${SPHEREDIR})
file (COPY round_lad.stl DESTINATION ${SPHEREDIR})
