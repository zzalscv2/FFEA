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

import sys


def readAndStore(iFile):
  H = []
  dH = {}
  with open(iFile, 'r') as sta:
    while (sta.readline() != "Measurements:\n"):
      continue

    blob = ""
    readblob = False
    for i in sta.readline().split():
      if i == "|":
        readblob = True
        continue
      if readblob:
        blob = i
        readblob = False
        continue
      H.append([blob+i])
    for line in sta:
      cnt = 0
      if line.count("RESTART"): continue
      for i in line.split():
        H[cnt].append(float(i))
        cnt += 1

  for i in range(len(H)):
    dH[H[i][0]] = i
    H[i].pop(0)

  return H, dH


iFile = ["sphere_63_120_two_measurement.out"]
H, dH = readAndStore(iFile[0])

## compare the PreCompEnergy of the first step to a value that has been checked:
err = 0
# if ( abs (H[dH["PreCompEnergy"]][-1] - 1.0075630000000001e-18) > 1e-32 ): 
if ( abs (H[dH["PreCompEnergy"]][-1] - 6.87781100000000001e-17) > 1e-32 ): 
  print( "PreComputed energy should be %e, but was found to be %e" % (1.0075630000000001e-18, (H[dH["PreCompEnergy"]][-1])))
  err = 1

sys.exit(err)
