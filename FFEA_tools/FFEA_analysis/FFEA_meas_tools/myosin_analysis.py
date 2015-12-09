# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 16:37:06 2015

@author: py12rw
"""

import argparse
import myosin_analysis_lib
import json
import numpy as np

def convert_string_list(string_list):
    """
    Converts a list of strings (e.g. ["3", "5", "6"]) to a list of integers.
    In: list of strings
    Out: list of integers
    """
    int_list = []
    for string in string_list:
        int_list.append(int(string))
    return int_list
    
def convert_array_in_dict_to_list(dict):
    """
    This turns a dictionary of arrays into a dictionary of lists. Useful as
    the json library cannot work with arrays.
    In: a dictionary containing at least one array
    Out: a dictionary containing no arrays
    """
    for key, value in dict.iteritems():
        if type(value) == type(np.zeros([1])):
            dict[key] = value.tolist()
    return dict

parser = argparse.ArgumentParser(description="CLI for analysing Myosin-7 trajectory files")

parser.add_argument("FFEA_filename", action="store", help="Path to input .ffea file")
parser.add_argument("-head_pin_file", action="store", help="Path to pinfile containing head")
parser.add_argument("-tail_pin_file", action="store", help="Path to pinfile containing tail")
parser.add_argument("-head_point_index", action="store", help="Index of a node in the center of the head")
parser.add_argument("-middle_point_index", action="store", help="Index of a node in the center of the molecule (ideally in the bendiest possible part)")
parser.add_argument("-tail_point_index", action="store", help="Index of a node in the center of the tail")
parser.add_argument("-twist", action="store_true", help="Whether to perform a twist test")
parser.add_argument("-head_line", action="store", nargs="+", help="Indices of points forming a line through the molecule\'s head, orthogonal to the principal axis (e.g. 3 34 51 57).")
parser.add_argument("-tail_line", action="store", nargs="+", help="Indices of points forming a line through the molecule\'s tail, orthogonal to the principal axis.")
parser.add_argument("-dist", action="store_true", help="Perform an end-to-end distance test (note: if this is disabled, the twist test will give results in terms of rads instead of rads/m)")
parser.add_argument("-angle", action="store_true", help="Perform an angular distribution test")
parser.add_argument("-plot", action="store_true", help="Save the results of any tests that have been performed to PDFs in the current working directory")
parser.add_argument("-name", action="store", help="Name of the molecule being tested (used in plots)")

args = parser.parse_args()

args.head_line = convert_string_list(args.head_line)
args.tail_line = convert_string_list(args.tail_line)

results = myosin_analysis_lib.run_sequential_tests(args.FFEA_filename,
                                                   args.head_pin_file,
                                                   args.tail_pin_file,
                                                   int(args.head_point_index),
                                                   int(args.tail_point_index),
                                                   int(args.middle_point_index),
                                                   args.head_line,
                                                   args.tail_line,
                                                   args.twist,
                                                   args.dist,
                                                   args.angle,
                                                   args.plot,
                                                   args.name)
                                   
results = convert_array_in_dict_to_list(results)

out_filename = args.FFEA_filename.split('.')[0]+"_myosin_analysis.json"

with open(out_filename, 'w') as outfile:
    json.dump(results, outfile)