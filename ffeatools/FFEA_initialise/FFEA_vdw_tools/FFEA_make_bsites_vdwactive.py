import sys, os
import FFEA_vdw, FFEA_binding_sites
import argparse as _argparse

parser = _argparse.ArgumentParser(description="FFEA Thin Trajectory")
parser.add_argument("bsitefname", action="store", help="INPUT .bsites")
parser.add_argument("vdwfname", action="store", help="INPUT .vdw)")
parser.add_argument("vdwoutfname", action="store", help="OUTPUT .vdw")
parser.add_argument("vdw_type", action="store", help="vdw type")

def make_bsites_vdwactive(bsite, vdwin, vdw_type):
    # Set all bsites as vdw active
    for s in bsites.bsites:
        for f in s.face_index:
            vdw.vdw_index[f] = vdw_type
    
    return vdw

if sys.stdin.isatty():
    args = parser.parse_args()
    # Get args and build objects

    # Build objects
    bsites = FFEA_binding_sites.FFEA_binding_sites(args.bsitefname)
    vdw = FFEA_vdw.FFEA_vdw(args.vdwinfname)

    # Write new file
    out_vdw = make_bsites_vdwactive(bsites, vdw, args.vdw_type)
    out_vdw.write_to_file(args.vdwoutfname)