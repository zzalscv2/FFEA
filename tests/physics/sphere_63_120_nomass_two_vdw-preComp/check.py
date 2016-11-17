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
if ( abs (H[dH["PreCompEnergy"]][-1] - 1.0075630000000001e-18) > 1e-32 ): 
  print( "PreComputed energy should be %e, but was found to be %e" % (1.0075630000000001e-18, (H[dH["PreCompEnergy"]][-1])))
  err = 1

sys.exit(err)
