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

from setuptools import setup, find_packages#, Extension
import subprocess
import sys as _sys

if len(_sys.argv) == 1:
    print("""
        This is the installer for the FFEA python package, but not FFEA itself! To install FFEA, do the following:
        
        cd ..
        mkdir FFEA_Build
        cd FFEA_Build
        cmake ../ffea
        make
        make install
        
        Alternatively, to install the FFEA python package, run 'python setup.py install'
        
        For more detailed information, please visit ffea.readthedocs.io
        """)
    _sys.exit(1)

def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) #dragons
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout

subprocess_cmd('mkdir ffeatools_build')
# there is a better way of compiling these, but doing so would disturb the all-inclusive build process
# todo: make these tools work as python extensions instead of being standalone
"""
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_analysis/Euler_characteristic_analysis; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/EM_density_map_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_node_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_mapping_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_volume_tools/make_cuboid_mesh; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_align; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')        
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')               
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/surface_coarse_grainer; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
"""

setup(name='ffeatools',
      version='1.0',
      description='FFEA input file generation and analysis tools',
      url='http://ffea.bitbucket.com',
      author='FFEA Team',
      author_email='???',
      license='GPLv3',
      packages=find_packages(),
      install_requires=[
          'numpy',
#          'pymol',
          'matplotlib',
          'argparse',
      ],
      zip_safe=False,
#      ext_modules = [
#          Extension('c_extension', sources = ['ffeatools/cpp_extension.cpp']),
#          Extension( )
#          ],
    )
# Todo: start converting C++ into boost.python stuff!
