
Tagging the binaries {#ffeaTags}
=================================================

This package comes with a mechanism to automatically tag the ffea executable
 with its exact version of the software, linking simulations to specific versions 
 of the code, with the aim to help users and developers in tracking bugs.

If at configure time "FFEA_DEVEL" i set to TRUE in CMake, 
 just the FFEA version (0.97, for example) will be printed at the begining of 
 the execution of "ffea". Otherwise, the commit tag, date, and branch 
 will be printed too. Details on how this is configured follow. Keep reading if 
 you are willing to modify the code.

At configure time, CMake uses "git" (the executable) to configure the sources
 so that the commit tag, date and branch of the version you're compiling are 
 printed at the beginning of the execution of "ffea". 

If "git" (the executable) is not found, or unusable (because at configure time 
 you don't have the repo but just the sources), the previous commit info
 (tag, date and branch) will be used. These details are written in a file
 (FFEA_version.h.in) that self-updates every time you do a commit. 

In order to do this automatic update of FFEA_version.h.in at commit time,
 "git" uses a "hook" (keep reading). There is a "pre-commit" hook in the FFEA root 
 folder, and in order to use it, every developer should create the following soft link:

     cd .git/hooks
     ln -s ../../pre-commit . 

The "pre-commit" hook is a simple bash script that uses "git" and "sed" to 
 configure FFEA_version.h.in, and add it to the commit.

