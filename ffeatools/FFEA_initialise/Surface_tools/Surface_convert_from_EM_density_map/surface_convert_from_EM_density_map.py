import wrap

def surface_convert_from_EM_density_map(*argv):
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./emdb_map_to_ffea", argv)
    return