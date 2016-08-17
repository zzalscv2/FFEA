import wrap

def build_euler_characteristic_table(*argv):
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./build_euler_characteristic_table", argv)
    return
