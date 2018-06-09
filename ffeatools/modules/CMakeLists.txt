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

install (FILES __init__.py wrap.py FFEA_meta.py
         FFEA_trajectory.py FFEA_measurement.py FFEA_frame.py 
         FFEA_binding_sites.py FFEA_kinetic_map.py FFEA_kinetic_rates.py
         FFEA_kinetic_states.py FFEA_material.py FFEA_node.py
         FFEA_pdb.py FFEA_pin.py FFEA_script.py FFEA_springs.py
         FFEA_stokes.py FFEA_surface.py FFEA_topology.py FFEA_turbotrajectory.py
         FFEA_vdw.py FFEA_universe.py FFEA_lj.py FFEA_exceptions.py
         FFEA_beads.py FFEA_ctforces.py FFEA_rod.py FFEA_skeleton.py
         DESTINATION "${PYTHONSTUFF}/modules")

