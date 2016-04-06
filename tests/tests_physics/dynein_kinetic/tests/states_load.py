import sys, os
import FFEA_kinetic_states

if len(sys.argv) != 2:
	sys.exit("Usage: python states_load.py [INPUT FFEA kinetic states (.states)]")

fname = sys.argv[1]

states = FFEA_kinetic_states.FFEA_kinetic_states(fname)
states.print_details()
