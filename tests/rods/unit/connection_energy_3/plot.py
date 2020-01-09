#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:00:00 2019

@author: rob
""" 

def read(log_file):
    results_dict = {0:{}, 1: {}, 2: {}, 3:{} }
    for key, value in results_dict.iteritems():
        results_dict[key] = { "+x": {}, "-x":{}, "+y":{}, "-y":{}, "+z":{}, "-z":{} }
        
    with open(log_file) as log:
        for line in log.readlines():
            if "ENERGYPLOT" in line:
                line_list = line.split(" ")
                displacement = float(line_list[2])
                node_index = int(line_list[4])
                results_dict[node_index]["+x"][displacement] = float(line_list[6])
                results_dict[node_index]["+y"][displacement] = float(line_list[7])
                results_dict[node_index]["+z"][displacement] = float(line_list[8])
                results_dict[node_index]["-x"][displacement] = float(line_list[9])
                results_dict[node_index]["-y"][displacement] = float(line_list[10])
                results_dict[node_index]["-z"][displacement] = float(line_list[11])
    return results_dict

def plot_one(plusdict, minusdict, node_index, axis, filename):
    x = []
    y = []
    minuskeys = minusdict.keys()
    minuskeys.reverse()
    minusvalues = minusdict.values()
    minusvalues.reverse()
    for key in minuskeys:
        x.append(key*-1)
    for key in plusdict.keys():
        x.append(key)
    for value in minusvalues:
        y.append(value)
    for value in plusdict.values():
        y.append(value)
        
    plt.plot(x, y)
    plt.xlabel("Displacement")
    plt.ylabel("Energy")
    plt.title("Node "+str(node_index)+" energies in the "+axis+" axis")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.clf()
    plt.close()
    return x, y

def plot_all(results_dict):
    for index in [0,1,2,3]:
        for axis in ["x", "y", "z"]:
            plot_one(results_dict[index]["+"+axis], results_dict[index]["-"+axis], index, axis, "node_"+str(index)+"_"+axis+"_axis.pdf")
    return

def auto():
    plot_all(read("/home/rob/software/FFEA_Build/tests/rods/unit/connection_energy_2/log.txt"))