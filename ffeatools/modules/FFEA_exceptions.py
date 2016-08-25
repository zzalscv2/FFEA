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
