import numpy as np
import sys

def readAndAverage(iFile, ini, end, th, Fields):
  H = []
  A = {}
  with open(iFile, 'r') as sta:
    while (sta.readline() != "Measurements:\n"):
      continue
    
    for i in sta.readline().split():
      H.append([i.strip()])
    for line in sta:
      cnt = 0
      if line.count("RESTART"): continue
      for i in line.split():
        H[cnt].append(float(i))
        cnt += 1

  if (ini == 0): ini = 1
  for i in H:
    if Fields.count(i[0]) == 0:
      continue
    print i[0], abs(np.mean(i[ini:end])/th - 1), np.std(i[ini:end]), len(i[ini:end])
    A[i[0]] = [abs(np.mean(i[ini:end])/th - 1), np.std(i[ini:end])]
  return A

nodes = 63
KbT = 4.11e-21
E = KbT*(3*nodes - 6)/2
print "Equipartition th: ", E
Tol = {"KineticEnergy": 0.09, "StrainEnergy":0.03}
ini = 40
end = -1

iFile = ["sphere_63_120_mass_measurement.out"]
         

err = 0
for f in iFile:
  print f
  A = readAndAverage(f, ini, end, E, Tol.keys())
  for s in Tol.keys():
    print A[s]
    if A[s][0] < Tol[s]: 
      print s, ": correct ", A[s][0], " < ", Tol[s]
    else: 
      print s, ": failed ", A[s][0], " > ", Tol[s]
      err = 1

  print "\n\n"


sys.exit(err)
