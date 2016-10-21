from setuptools import setup, find_packages#, Extension
import subprocess

def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) #dragons
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout

subprocess_cmd('mkdir ffeatools_build')
# there is a better way of compiling these, but doing so would disturb the all-inclusive build process
# todo: make these tools work as python extensions instead of being standalone
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_analysis/Euler_characteristic_analysis; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/EM_density_map_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_node_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_mapping_tools; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_volume_tools/make_cuboid_mesh; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_align; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')        
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')               
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/surface_coarse_grainer; make; make install; rm Makefile CMakeCache.txt cmake_install.cmake')
               
setup(name='ffeatools',
      version='1.0',
      description='FFEA file generation and analysis tools',
      url='http://ffea.bitbucket.com',
      author='FFEA Team',
      author_email='???',
      license='???',
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