from setuptools import setup#, Extension
import subprocess

def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) #dragons
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout

# there is a better way of compiling these, but doing so would disturb the all-inclusive build process
# todo: make these tools work as python extensions instead of being standalone
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_analysis/Euler_characteristic_analysis; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/EM_density_map_tools; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_node_tools; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_mapping_tools; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/FFEA_volume_tools/make_cuboid_mesh; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_align; make; make install; bash clean.sh')        
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/PDB_tools/PDB_convert_from_FFEA_trajectory; make; make install; bash clean.sh')
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map; make; make install; bash clean.sh')               
subprocess_cmd('cd ffeatools_build; cmake ../ffeatools/FFEA_initialise/Surface_tools/surface_coarse_grainer; make; make install; bash clean.sh')
               
setup(name='ffeatools',
      version='0.9',
      description='FFEA file generation and analysis tools',
      url='http://ffea.bitbucket.com',
      author='FFEA Team',
      author_email='???',
      license='???',
      packages=['ffeatools'],
      install_requires=[
          'numpy',
#          'pymol',
          'matplotlib',
          'argparse',
#          'math'
      ],
      zip_safe=False,
#      ext_modules = [
#          Extension('c_extension', sources = ['ffeatools/cpp_extension.cpp']),
#          Extension( )
#          ],
    )
# When we start converting the C++ into boost.python stuff, I will uncomment this.