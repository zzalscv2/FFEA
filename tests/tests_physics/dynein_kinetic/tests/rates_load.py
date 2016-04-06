import sys, os
import FFEA_kinetic_rates

if len(sys.argv) != 2:
	sys.exit("Usage: python rates_load.py [INPUT FFEA kinetic rates (.rates)]")

fname = sys.argv[1]

rates = FFEA_kinetic_rates.FFEA_kinetic_rates(fname)
rates.print_details()
