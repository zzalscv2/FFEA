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

class FFEAFormatError(Exception):

	def __init__(self, lin="Unknown", lstr=""):

		self.lin = str(lin)
		self.lstr = lstr

class FFEAIOError(IOError):

	def __init__(self, fname="", fext=[""]):

		self.fname = str(fname)
		self.fext = fext

def print_error():
	from sys import stdout
	stdout.write("\n\033[91mERROR:\033[0m ")
