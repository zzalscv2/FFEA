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

import sys, os

if len(sys.argv) < 5:
	sys.exit("Usage: python " + sys.argv[0] + " [spring_file_fname] [spring_constant] [tolerance_length] [List of pairs of nodes]")

spring_file_fname = sys.argv[1]
k = float(sys.argv[2])
l = float(sys.argv[3])
node_pairs = []
for i in range(4, len(sys.argv), 2):
	node_pairs.append([sys.argv[i], sys.argv[i + 1]])

fout = open(spring_file_fname, "w")
fout.write("ffea spring file\n")
fout.write("num_springs %d\n" % (len(node_pairs)))
fout.write("springs:\n")
for node_pair in node_pairs:
	fout.write("0 0 %d 1 0 %d %e %e\n" % (int(node_pair[0]), int(node_pair[1]), k, l))

fout.close()
