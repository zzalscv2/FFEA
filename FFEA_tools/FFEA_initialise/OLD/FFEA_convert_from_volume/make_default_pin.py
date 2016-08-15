import os, sys

if len(sys.argv) != 2:
	sys.exit("make_default_pin.py [OUTPUT PIN FNAME]")

print "Running: make_default_pin.py"

pin_txt = "walrus pinned nodes file\nnum_pinned_nodes 0\npinned nodes:\n"

out_fname = sys.argv[1]
with open(out_fname, "w") as outfile:
	outfile.write(pin_txt)

print "Done. --> make_default_pin.py"
