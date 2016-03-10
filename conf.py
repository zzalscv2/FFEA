import subprocess, os, re

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

# read_the_docs_build = True

if read_the_docs_build:

    sta = open("Doxygen.in",'r')
    STA = sta.readlines()
    sta.close()
   

    sta = open("Doxyfile",'w')
    for i in STA: 
      txt = re.sub("@CMAKE_CURRENT_SOURCE_DIR@",".",i)
      txt = re.sub("@PACKAGE_VERSION@","0.9",txt)
      l = txt.split()
      if len(l) > 1:
        if l[0] == "OUTPUT_DIRECTORY":
          txt = "OUTPUT_DIRECTORY		= .\n"
      sta.write(txt)
    sta.close()
    
    subprocess.call('cd ../doxygen; doxygen', shell=True)
