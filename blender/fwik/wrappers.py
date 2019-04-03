import bpy
import mathutils
import math

from . import fwik

class FWIKBone:
    def __init__(self, armature, pose_bone):
        self.armature = armature
        self.pose_bone = pose_bone

    def get_alpha_position(self, alpha):
        return self.get_tail_position() * alpha + self.get_head_position() * (1 - alpha)

    ############################
    #### Simulation Helpers ####
    ############################

    def get_key(self):
        return self.pose_bone.name

    def get_parent(self):
        if self.pose_bone.parent == None:
            return None
        else:
            bone = self.armature.pose.bones[self.pose_bone.parent.name]
            return FWIKBone(self.armature, bone)

    def get_average_radius(self):
        return 0.5

    def get_child(self):
        if self.pose_bone.child == None:
            return None
        else:
            bone = self.armature.pose.bones[self.pose_bone.child.name]
            return FWIKBone(self.armature, bone)

    def get_tail_position(self):
        return self.armature.matrix_world * self.pose_bone.tail

    def get_head_position(self):
        return self.armature.matrix_world * self.pose_bone.head

    def get_length(self):
        return (self.get_head_position() - self.get_tail_position()).length

    def get_axial_rotation(self):
        mode = self.pose_bone.rotation_mode
        if mode == 'QUATERNION':
            return self.pose_bone.rotation_quaternion.copy()
        elif mode == 'AXIS_ANGLE':
            return mathutils.Quaternion(self.pose_bone.rotation_axis_angle)
        else:
            return self.pose_bone.rotation_euler.to_quaternion()

    def set_axial_rotation(self, lr):
        mode = self.pose_bone.rotation_mode
        if mode == 'QUATERNION':
            self.pose_bone.rotation_quaternion = lr
        elif mode == 'AXIS_ANGLE':
            self.pose_bone.rotation_axis_angle = lr.to_axis_angle()
        else:
            self.pose_bone.rotation_euler = lr.to_euler(mode)

    def get_world_rotation(self):
        return (self.armature.matrix_world * self.pose_bone.matrix).to_quaternion()

    def get_bone_axes_to_world(self):
        x = self.pose_bone.x_axis
        y = self.pose_bone.y_axis
        z = self.pose_bone.z_axis

        bone_axes = mathutils.Matrix([
            [ x.x, y.x, z.x ],
            [ x.y, y.y, z.y ],
            [ x.z, y.z, z.z ]
        ])

        return self.armature.matrix_world.to_3x3() * bone_axes

    def get_mass(self):
        return self.pose_bone['fwik_mass']

    def get_angular_velocity(self):
        return mathutils.Vector(self.pose_bone['fwik_angular_velocity'])

    def set_angular_velocity(self, av):
        self.pose_bone['fwik_angular_velocity'] = av

    def get_linear_velocity(self):
        return mathutils.Vector(self.pose_bone['fwik_linear_velocity'])

    def set_linear_velocity(self, lv):
        self.pose_bone['fwik_linear_velocity'] = lv

    def use_ik_limit_x(self):
        return self.pose_bone.use_ik_limit_x

    def get_min_rot_x(self):
        return self.pose_bone.ik_min_x

    def get_max_rot_x(self):
        return self.pose_bone.ik_max_x

    def use_ik_limit_y(self):
        return self.pose_bone.use_ik_limit_y

    def get_min_rot_y(self):
        return self.pose_bone.ik_min_y

    def get_max_rot_y(self):
        return self.pose_bone.ik_max_y

    def use_ik_limit_z(self):
        return self.pose_bone.use_ik_limit_z

    def get_min_rot_z(self):
        return self.pose_bone.ik_min_z

    def get_max_rot_z(self):
        return self.pose_bone.ik_max_z

    def translate(self, translation):
        self.pose_bone.location += translation

class FWIKControlPoint:
    def __init__(self, object):
        assert(FWIKControlPoint.isFWIKControlPoint(object))
        self.object = object

    ###########################
    #### Interface Methods ####
    ###########################

    @property
    def armature(self):
        return self.object['fwik_armature']

    @property
    def bone_name(self):
        return self.object['fwik_bone_name']

    @property
    def active(self):
        return self.object['fwik_active']

    @property
    def bone_alpha(self):
        return self.object['fwik_bone_alpha']

    ############################
    #### Simulation Helpers ####
    ############################

    def get_position(self):
        return self.object.location

    def get_bone(self):
        return FWIKBone(self.armature, self.armature.pose.bones[self.bone_name])

    def get_attachment_position(self):
        return self.get_bone().get_alpha_position(self.bone_alpha)

    def get_spring_constant(self):
        return self.object['fwik_spring_constant']

    ##################################
    #### Static Interface Methods ####
    ##################################

    def create(parent_rig, armature, bone, active=False, bone_alpha=1.0, spring_constant=10.0):
        # Start with an empty object
        control_point = bpy.data.objects.new("FWIK Control Point", None)
        control_point.empty_draw_size = 0.1

        # Keep track of the fact that it is a fwik rig with a read-only property
        control_point['fwik_is_control_point'] = True

        # Other properties
        control_point['fwik_parent_rig'] = parent_rig
        control_point['fwik_armature'] = armature
        control_point['fwik_bone_name'] = bone.name

        control_point['fwik_active'] = active
        control_point['fwik_bone_alpha'] = bone_alpha
        control_point['fwik_spring_constant'] = spring_constant

        fcp = FWIKControlPoint(control_point)
        fcp.object.location = fcp.get_bone().get_alpha_position(bone_alpha)

        return fcp

    def isFWIKControlPoint(object):
        return object != None and getattr(object, 'fwik_is_control_point', False)

class FWIKRig:

    def __init__(self, object):
        assert(FWIKRig.isFWIKRig(object))
        self.object = object
    
    ###########################
    #### Interface Methods ####
    ###########################

    @property
    def armature(self):
        return self.object['fwik_armature']

    @property
    def iterations(self):
        return self.object['fwik_iterations']

    @property
    def damping(self):
        return self.object['fwik_damping']

    ############################
    #### Simulation Methods ####
    ############################

    def get_control_points(self):
        return [FWIKControlPoint(child) for child in self.object.children if FWIKControlPoint.isFWIKControlPoint(child)]

    def get_bones(self):
        return [FWIKBone(self.armature, bone) for bone in self.armature.pose.bones]

    def simulate(self, iterations, dt):
        fwik.Simulator(rig=self, iterations=iterations, time_step=dt).run()

    ##################################
    #### Static Interface Methods ####
    ##################################

    def create(armature, iterations=50, damping=3.0):
        # Start with an empty object
        rig = bpy.data.objects.new("FWIK Rig", None)

        # Make it "invisible"
        rig.empty_draw_size = 0

        # Keep track of the fact that it is a fwik rig with a read-only property
        rig['fwik_is_rig'] = True

        # Other properties
        rig['fwik_armature'] = armature
        rig['fwik_iterations'] = iterations
        rig['fwik_damping'] = damping

        return FWIKRig(rig)

    def isFWIKRig(object):
        return object != None and getattr(object, 'fwik_is_rig', False)

class MakeFWIKControlPoint(bpy.types.Operator):
    """Creates a control point for a FWIK simulation"""
    bl_idname = "object.make_fwik_control_point"
    bl_label = "Make FWIK Control Point"
    bl_options = {'REGISTER', 'UNDO'} 

    is_active = bpy.props.BoolProperty(name="Active", default=True)
    bone_alpha = bpy.props.FloatProperty(name="Bone Alpha",\
        description="How far along the bone should the control point be attached to",\
        min=0.0, max=1.0, default=1.0)
    spring_constant = bpy.props.FloatProperty(name="Spring Constant",\
        description="The spring constant of the 'spring' connected to the bone. Higher means more influence.",\
        min=0.0, default=10.0)

    @classmethod
    def poll(cls, context):
        return context.mode == 'POSE'

    def execute(self, context):
        scene = context.scene

        active = scene.objects.active

        # Do we have a rig selected?
        fwik_rig = None
        armature = None
        bone = None

        if FWIKRig.isFWIKRig(active):
            fwik_rig = active
            armature = fwik_rig['fwik_armature']

        if active != None and FWIKRig.isFWIKRig(active.parent):
            fwik_rig = active.parent
            armature = active

        if fwik_rig == None:
            self.report({'ERROR'}, 'Not in a FWIK Rig. Select an armature and use "Make FWIG Rig" first.')
            return {'CANCELLED'}

        if armature != None:
            bone = armature.data.bones.active

        #####################################
        #### Create Control Point object ####
        #####################################

        control_point = FWIKControlPoint.create(parent_rig=fwik_rig, armature=armature, bone=bone,\
            active=self.is_active, bone_alpha=self.bone_alpha, spring_constant=self.spring_constant)

        ######################
        #### Put in scene ####
        ######################

        if fwik_rig != None:
            control_point.object.parent = fwik_rig

        scene.objects.link(control_point.object)

        return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self)

    def draw(self, context):
        layout = self.layout
        layout.prop(self, 'is_active')
        layout.prop(self, 'bone_alpha')
        layout.prop(self, 'spring_constant')

class MakeFWIKRig(bpy.types.Operator):
    """Creates a rig for a FWIK simulation"""
    bl_idname = "object.make_fwik_rig"
    bl_label = "Make FWIK Rig"
    bl_options = {'REGISTER', 'UNDO'} 

    def execute(self, context):
        scene = context.scene

        active = scene.objects.active

        # Ensure we have an armature selected
        if active == None or active.type != "ARMATURE":
            self.report({'ERROR'}, 'No armature selected.')
            return {'CANCELLED'}

        # Make sure its not already in a FWIK Rig
        if active != None and FWIKRig.isFWIKRig(active.parent):
            self.report({'ERROR'}, 'Armature is already in a FWIK Rig.')
            return {'CANCELLED'}

        ###########################
        #### Create rig object ####
        ###########################

        rig = FWIKRig.create(active)

        ######################
        #### Put in scene ####
        ######################

        scene.objects.link(rig.object)
        active.parent = rig.object
    
        # Set custom props
        for bone in active.pose.bones:
            bone['fwik_mass'] = 1.0

        return {'FINISHED'}

class TestFWIK(bpy.types.Operator):
    """Runs a hardcoded test for FWIK"""
    bl_idname = "object.test_fwik"
    bl_label = "Test FWIK"

    _timer = None

    def execute(self, context):
        scene = context.scene
        active = scene.objects.active

        fwik_rig = None

        if FWIKRig.isFWIKRig(active):
            fwik_rig = active

        if active != None and FWIKRig.isFWIKRig(active.parent):
            fwik_rig = active.parent

        if fwik_rig != None:
            self.frig = FWIKRig(fwik_rig)

            for bone in self.frig.get_bones():
                bone.pose_bone['fwik_angular_velocity'] = mathutils.Vector((0,0,0))
                bone.pose_bone['fwik_linear_velocity'] = mathutils.Vector((0,0,0))

            self._timer = context.window_manager.event_timer_add(.1, context.window)
            context.window_manager.modal_handler_add(self)
            print('Running animation...')
            return {'RUNNING_MODAL'}
        else:
            print('No FWIK rig selected.')
            return {'CANCELLED'}


    def modal(self, context, event):
        if event.type == 'ESC':
            return self.cancel(context)

        if event.type == 'TIMER':
            self.frig.simulate(iterations=1, dt=0.05)

        return {'PASS_THROUGH'}

    def cancel(self, context):
        context.window_manager.event_timer_remove(self._timer)
        print('Stopped')
        return {'CANCELLED'}

########################
#### Register Props ####
########################

bpy.types.Object.fwik_is_control_point = bpy.props.BoolProperty(name="Is FWIK Control Point",\
    default=False)
bpy.types.Object.fwik_parent_rig = bpy.props.PointerProperty(name="Parent Rig",\
    type=bpy.types.Object)
bpy.types.Object.fwik_armature = bpy.props.PointerProperty(name="Armature",\
    type=bpy.types.Object)
bpy.types.Object.fwik_bone_name = bpy.props.StringProperty(name="Bone")
bpy.types.Object.fwik_active = bpy.props.BoolProperty(name="Active")
bpy.types.Object.fwik_bone_alpha= bpy.props.FloatProperty(name="Bone Alpha",\
    min=0.0, max=1.0)
bpy.types.Object.fwik_spring_constant = bpy.props.FloatProperty(name="Spring Constant",\
    min=0.0)
bpy.types.Object.fwik_is_rig = bpy.props.BoolProperty(name="Is FWIK Rig",\
    default=False)
bpy.types.Object.fwik_iterations = bpy.props.IntProperty(name="Iterations",\
    default=50, min=1)
bpy.types.Object.fwik_damping  = bpy.props.FloatProperty(name="Damping",\
    min=0.0)