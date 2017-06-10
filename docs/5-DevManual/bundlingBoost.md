
How to bundle Boost {#bundleBoost}
=================================================

Currently, the code is shipped with a subset of Boost v. 1.63. In the future, 
 somebody may think of upgrading to a new version, or of using more Boost modules.
 That is how I took the current subset: 

       tar -jxf boost_1_63_0.tar.bz2
       cd boost_1_63_0
       ./bootstrap.sh
       ./bjam tools/bcp
       ls ./tools/bcp/
       mkdir $HOME/install/boost/1.63_bundle
       ./dist/bin/bcp tools/build boost/program_options.hpp boost/algorithm/string.hpp boost/filesystem.hpp boost/lexical_cast.hpp $HOME/install/boost/1.63_bundle
       ./dist/bin/bcp --report tools/build boost/program_options.hpp boost/algorithm/string.hpp boost/filesystem.hpp boost/lexical_cast.hpp $HOME/install/boost/1.63_bundle.html

 
 

