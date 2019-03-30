bl_info = {
    "name": "FWIK",
    "description": "A utility to calculate inverse kinematics by forward physical simulation.",
    "author": "koreanwglasses",
    "version": (0, 0, 0),
    "warning": "This addon is still in development.",
    "category": "Animation" 
}


import bpy


# load and reload submodules
##################################

import importlib
from . import developer_utils
importlib.reload(developer_utils)
modules = developer_utils.setup_addon_modules(__path__, __name__, "bpy" in locals())



# register
##################################

import traceback

def register():
    try: bpy.utils.register_module(__name__)
    except: traceback.print_exc()

    print("Registered {} with {} modules".format(bl_info["name"], len(modules)))

def unregister():
    try: bpy.utils.unregister_module(__name__)
    except: traceback.print_exc()

    print("Unregistered {}".format(bl_info["name"]))
