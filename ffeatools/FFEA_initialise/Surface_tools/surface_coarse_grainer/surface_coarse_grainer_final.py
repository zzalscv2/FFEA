import wrap

def surface_coarse_grainer(*argv):
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./surface_coarse_grainer_final", argv)
    return