#!/bin/csh -fx
#
# test of bin/pcazip utilities
#
\rm -f testm.pcz zippedm.pcz short1.pcz short2.pcz albumm.pcz long.mdcrd long.dcd long.pcz
../bin/pcazip -help
../bin/pcazip -i test.x -n 758 -o testm.pcz -v -mask mask.pdb
../bin/pcazip -i short.x -n 758 -o short1.pcz -v
../bin/pcazip -i short.x -n 758 -o short2.pcz -nofast -v
../bin/pcazip -a test.alb -n 758 -o albumm.pcz -v -mask mask.pdb -v
../bin/pcaunzip -i albumm.pcz -o long.mdcrd
../bin/pcaunzip -i albumm.pcz -o long.dcd -format charmm
head -2 long.mdcrd
../bin/pcaunzip -i albumm.pcz -iv1 2 -iv2 4 |  head -2
../bin/pcazip -i long.dcd -n 22 -o long.pcz -q 99
../bin/pczdump -help
../bin/pczdump -info -i testm.pcz
../bin/pczdump -avg -i testm.pcz
../bin/pczdump -avg -pdb mask.pdb -i testm.pcz
../bin/pczdump -evals -i testm.pcz
../bin/pczdump -evec 1 -i testm.pcz
../bin/pczdump -rms 1 -i testm.pcz
../bin/pczdump -fluc 1 -i testm.pcz
../bin/pczdump -anim 2 -i testm.pcz
../bin/pczdump -anim 2 -pdb mask.pdb -i testm.pcz
../bin/pczdump -proj 1 -i testm.pcz
../bin/pczdump -coll -i testm.pcz
../bin/pczdump -maha 10 -i testm.pcz
../bin/pczcomp -x testm.pcz -y albumm.pcz
../bin/pczcomp -x long.pcz -y albumm.pcz
