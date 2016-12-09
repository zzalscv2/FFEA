import wrap
import __builtin__

def surface_convert_from_EM_density_map(*argv):
    argv = wrap.sanitize_tuple(argv)
    
    # Log
    try:
        __builtin__.meta_info.log["Surface convert from EM density map"] = argv
    except:
        print("Warning! failed to write log.")
        pass    
    
    wrap.wrap_process("./emdb_map_to_ffea", argv)
    return