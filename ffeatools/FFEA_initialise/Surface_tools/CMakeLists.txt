
add_subdirectory("surface_coarse_grainer")
add_subdirectory("Surface_convert_from_EM_density_map")

install(FILES __init__.py
        "FFEA_surf_and_node_convert_from_obj.py" "Gts_convert_from_netgen_surf.py"
        "Gts_convert_from_obj.py" "Netgen_surf_convert_from_obj.py"
        "Netgen_surf_convert_from_off.py" "Obj_convert_from_netgen_surf.py"
        "Off_convert_from_obj.py" "scale_netgen_surf.py"
        DESTINATION "${PYTHONSTUFF}/FFEA_initialise/Surface_tools")

