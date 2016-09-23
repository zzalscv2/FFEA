import wrap
import __builtin__

def surface_coarse_grainer(*argv):

    # Log
    try:
        surface_coarse_grainer = {}
        surface_coarse_grainer["Input .surf file"] = argv[0]
        surface_coarse_grainer["Output file"] = argv[1]
        surface_coarse_grainer["Coarseness level"] = argv[2]
        surface_coarse_grainer["Coarseness level"] = argv[3]
        surface_coarse_grainer["Conserve volume"] = argv[4]
        surface_coarse_grainer["Find smallest edge"] = argv[5]

        if len(argv) > 6:
            surface_coarse_grainer["Xmin"] = argv[6]
            surface_coarse_grainer["Xmax"] = argv[7]
            surface_coarse_grainer["Zmin"] = argv[8]
            surface_coarse_grainer["Zmax"] = argv[9]
            for i in range(100):
                if "Surface Coarse Grainer "+str(i) in __builtin__.meta_info.log:
                    pass
                else:
                    __builtin__.meta_info.log["Surface Coarse Grainer "+str(i+1)] = surface_coarse_grainer
                    break
        else:
            __builtin__.meta_info.log["Surface Coarse Grainer"] = surface_coarse_grainer
                          
    except:
        print("Warning! failed to write log.")
        pass        
    
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./surface_coarse_grainer_final", argv)
    return