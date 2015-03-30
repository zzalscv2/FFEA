import sys, os

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .pcz file] [INPUT reference topology (_frame0.pdb)] [num_animations]")

inpcz = sys.argv[1]
inref = sys.argv[2]
num_anim = int(sys.argv[3])
anim_basename = inpcz.split(".")[0]
for i in range(num_anim):
	os.system("pyPczdump --input " + inpcz + " --anim " + str(i) + " --pdb " + inref + " -o " + anim_basename + "_anim" + str(i) + ".pdb")
