#!/bin/bash
# This script will put the pyPCAZIP suite through its paces, exemplifying
# a lot of the most common ways it can be used.
echo " "
echo "pyPcazip: usage"
echo " "
pyPcazip -h
echo " "
echo "test 1: very basic"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o 2ozq.pcz --selection "name CA" -vvv -p justCA.pdb
echo " "
echo "test 2: with slicing of the trajectory file"
echo " "
pyPcazip --topology 2ozq.pdb -i '2ozq.dcd(1:10)' -o 2ozq.pcz --selection "name CA" -vvv
echo " "
echo "test 3: with two trajectory files"
echo " "
pyPcazip --topology 2ozq.pdb -i '2ozq.dcd(1:10)' '2ozq.xtc(::3)' -o 2ozq.pcz --selection "name CA" -vvv
echo " "
echo "test 4: with a mask file to select atoms to include"
echo " "
pyPcazip --topology 2ozq.pdb -i '2ozq.dcd(1:10)' -o 2ozq.pcz --mask pocket.pdb -vvv
echo " "
echo "test 5:  with an album"
echo " "
pyPcazip --topology 2ozq.pdb -a 2ozq.alb -o 2ozq.pcz --mask pocket.pdb -vvv
echo " "
echo "pyPcaunzip - usage:"
echo " "
echo "test 6:  using PZC7 version"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o 2ozq_PCZ7.pcz --mask pocket.pdb -f PCZ7 -vvv
echo " "

echo "pyPcaunzip - usage:"
echo " "
pyPcaunzip -h
echo " "
echo "test 7: Unzipping..."
echo " "
pyPcaunzip --topology pocket.pdb --compressed 2ozq.pcz -o 2ozq_proteinpocket.dcd -vvv
echo " "
echo "test 8: Unzipping PCZ7"
echo " "
pyPcaunzip --topology pocket.pdb --compressed 2ozq_PCZ7.pcz -o 2ozq_proteinpocket_unzipped.dcd -vvv
echo " "

echo "pyPczdump - usage:"
echo " "
pyPczdump -h
echo " "
echo "test 9: Print basic information from a compressed file:"
echo " "
pyPczdump --input 2ozq.pcz --info -vvv
echo " "
echo "test 10: Print the average structure out of a compressed trajectory file:"
echo " "
pyPczdump --input 2ozq.pcz --avg --pdb pocket.pdb -vvv
echo " "
echo "test 11: Print out the eigenvectors in a compressed file."
echo " "
pyPczdump --input 2ozq.pcz --evals -vvv
echo " "
echo "test 12: Print out a specific eigenvector."
echo " "
pyPczdump --input 2ozq.pcz --evec 1 -vvv
echo " "
echo "test 13: Print out the projections of a specific eigenvector."
echo " "
pyPczdump --input 2ozq.pcz --proj 1 -vvv
echo " "
echo "test 14: Print out the atomic fluctuations related to a specific eigenvector."
echo " "
pyPczdump --input 2ozq.pcz --fluc 1 -vvv
echo " "
echo "test 15: Produce an animation along eigenvector 1:"
echo " "
pyPczdump --input 2ozq.pcz --anim 1 --pdb pocket.pdb -vvv
echo " "
echo "test 16: Print out the rmsd of each frame from frame 1 (the second, as pyPCAZIP counts from 0)."
echo " "
pyPczdump --input 2ozq.pcz --rms 1 -vvv
echo " "
echo "test 17: Print out a collectivity metric for each eigenvector."
echo " "
pyPczdump --input 2ozq.pcz --coll -vvv
echo " "

echo "pyPczcomp - usage:"
echo " "
pyPczcomp -h
echo " "
echo "Basic example of pyPczcomp:"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.ncdf -o 2ozq_A.pcz --selection "name CA" -f PCZ4
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o 2ozq_B.pcz --selection "name CA" -f PCZ4
pyPcazip --topology 2ozq.pdb -i 2ozq.ncdf -o 2ozq_C.pcz --selection "name CA" -f PCZ6
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o 2ozq_D.pcz --selection "name CA" -f PCZ6
pyPcazip --topology 2ozq.pdb -i 2ozq.ncdf -o 2ozq_E.pcz --selection "name CA" -f PCZ7
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o 2ozq_F.pcz --selection "name CA" -f PCZ7
pyPczcomp -i 2ozq_A.pcz 2ozq_B.pcz
pyPczcomp -i 2ozq_A.pcz 2ozq_C.pcz
pyPczcomp -i 2ozq_A.pcz 2ozq_D.pcz
pyPczcomp -i 2ozq_A.pcz 2ozq_E.pcz
pyPczcomp -i 2ozq_A.pcz 2ozq_F.pcz

pyPczcomp -i 2ozq_B.pcz 2ozq_C.pcz
pyPczcomp -i 2ozq_B.pcz 2ozq_D.pcz
pyPczcomp -i 2ozq_B.pcz 2ozq_E.pcz
pyPczcomp -i 2ozq_B.pcz 2ozq_F.pcz

pyPczcomp -i 2ozq_C.pcz 2ozq_D.pcz
pyPczcomp -i 2ozq_C.pcz 2ozq_E.pcz
pyPczcomp -i 2ozq_C.pcz 2ozq_F.pcz

pyPczcomp -i 2ozq_D.pcz 2ozq_E.pcz
pyPczcomp -i 2ozq_D.pcz 2ozq_F.pcz

pyPczcomp -i 2ozq_E.pcz 2ozq_F.pcz
