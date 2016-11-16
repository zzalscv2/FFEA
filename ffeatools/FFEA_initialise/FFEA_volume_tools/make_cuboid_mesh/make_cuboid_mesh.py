import wrap

def make_cuboid_mesh(*argv):
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./make_cuboid_mesh", argv)
    return