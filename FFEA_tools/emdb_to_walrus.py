import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python emdb_to_walrus.py [EMDB ID NUMBER] [LEVEL]\n");

id = sys.argv[1]
level = sys.argv[2]
emdb_repo = "ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/"

emdfileloc = "EMD-" + id + "/map/"

emdb_file = "emd_" + id + ".map"

print "Getting emd file..."
os.system("wget " + emdb_repo + emdfileloc + emdb_file + ".gz\n")
print "...done"

print "Unzipping..."
os.system("gunzip " + emdb_file + ".gz\n")
print "...done"

print "Creating mesh from emdb map at level = " + level + "..."

print "...done"
