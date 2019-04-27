import bpy
from mathutils import *
import math
from itertools import combinations

class PhysBone:
    '''
    A wrapper for the bone to perform calculations and modify the actual bone
    '''
    def __init__(self, parent, bone):
        self.parent = parent
        self.bone = bone

        # Physics parameters
        self.new_local_axial_rotation = self.bone.get_axial_rotation()
        self.new_angular_velocity = self.bone.get_angular_velocity()
        self.new_linear_velocity = self.bone.get_linear_velocity()

        self.net_force = Vector((0,0,0))
        self.net_torque = Vector((0,0,0))

        # Precompute stuff
        self.mass = self.bone.get_mass()

        # Compute principle moments asusming bone is a cylinder
        m = self.mass
        r = self.bone.get_average_radius()
        h = self.bone.get_length()
        principle_moment_x = (1.0 / 12.0) * m * (3 * r**2 + h**2)
        principle_moment_y = (1.0 / 2.0) * m * r**2
        self.principle_moment = Matrix([[principle_moment_x, 0, 0 ], 
                                    [0, principle_moment_y, 0],
                                    [0, 0, principle_moment_x]])

        self.translation = Vector([0,0,0])

        self.memos = [None, None]
        self.integration_memos = [None, None, None, None, None]
        self.apply_memos = [None]


    def get_parent_bone(self):
        key = 0
        if self.memos[key] == None:
            self.memos[key] = self.parent.get_physbone(self.bone.get_parent())
        return self.memos[key]

    def get_child_bones(self):
        key = 1
        if self.memos[key] == None:
            self.memos[key] = [self.parent.get_physbone(bone) for bone in self.bone.get_children()]
        return self.memos[key]

    def get_head_position(self):
        key = 0
        if self.integration_memos[key] == None:
            self.integration_memos[key] = self.bone.get_head_position()
        return self.integration_memos[key]

    def get_tail_position(self):
        key = 1
        if self.integration_memos[key] == None:
            self.integration_memos[key] = self.bone.get_tail_position()
        return self.integration_memos[key]

    def get_center_position(self):
        key = 4
        if self.integration_memos[key] == None:
            self.integration_memos[key] = 0.5 * (self.get_head_position() + self.get_tail_position())
        return self.integration_memos[key]
    
    def get_world_rotation(self):
        return self.bone.get_world_rotation()

    def compute_head_velocity(self):
        parent_bone = self.get_parent_bone()
        if parent_bone == None:
            center_pos = self.get_center_position()
            head_pos = self.get_head_position()
            return self.new_linear_velocity + self.new_angular_velocity.cross(head_pos - center_pos)
        else:
            return parent_bone.compute_tail_velocity()

    def compute_tail_velocity(self):
        tail_pos = self.get_tail_position()
        head_pos = self.get_head_position()
        head_vel = self.compute_head_velocity()
        return head_vel + self.new_angular_velocity.cross(tail_pos - head_pos)
    
    def compute_angular_acceleration(self):
        key = 0
        if self.apply_memos[key] == None:
            Icm = self.compute_moment()
            self.apply_memos[key] = Icm.inverted() * (self.net_torque - self.new_angular_velocity.cross(Icm * self.new_angular_velocity))
        return self.apply_memos[key]

    def compute_linear_acceleration(self):
        return self.net_force / self.mass

    def compute_head_acceleration(self):
        linear_accel = self.compute_linear_acceleration()
        angular_accel = self.compute_angular_acceleration()
        head_pos = self.get_head_position()
        center_pos = self.get_center_position()
        return linear_accel + angular_accel.cross(head_pos - center_pos)

    def compute_tail_acceleration(self):
        linear_accel = self.compute_linear_acceleration()
        angular_accel = self.compute_angular_acceleration()
        tail_pos = self.get_tail_position()
        center_pos = self.get_center_position()
        return linear_accel + angular_accel.cross(tail_pos - center_pos)

    def compute_acceleration_at(self, point):
        linear_accel = self.compute_linear_acceleration()
        angular_accel = self.compute_angular_acceleration()
        center_pos = self.get_center_position()

        return linear_accel + angular_accel.cross(point - center_pos)

    def compute_force_acceleration_matrix(self, point):
        r = point - self.get_center_position()
        rx = Matrix([[0, -r.z, r.y],
                     [r.z, 0, -r.x],
                     [-r.y, r.x, 0]])

        return Matrix.Scale(1.0 / self.mass, 3) + rx.transposed() * self.compute_moment().inverted() * rx

    def compute_moment(self):
        key = 2
        if self.integration_memos[key] == None:
            P = self.get_world_to_axial()
            self.integration_memos[key] = P.inverted() * self.principle_moment * P
        return self.integration_memos[key]

    def compute_parent_angular_velocity(self):
        parent = self.get_parent_bone()
        if parent == None:
            return Vector([0, 0, 0])
        else:
            return parent.new_angular_velocity

    def reset(self):
        self.net_force = Vector((0,0,0))
        self.net_torque = Vector((0,0,0))

    def apply_force(self, force, position):
        self.net_force += force

        r = position - self.get_center_position()
        torque = r.cross(force)

        self.net_torque += torque

        self.apply_memos = [None]

    def apply_torque(self, torque):
        self.net_torque += torque
        self.apply_memos = [None]

    def apply_axial_torque(self, axial_torque):
        self.net_torque += self.get_world_to_axial().inverted() * axial_torque
        self.apply_memos = [None]

    def get_world_to_axial(self):
        key = 3
        if self.integration_memos[key] == None:
            bone_axes = self.bone.get_bone_axes_to_world()
            axial = self.bone.get_axial_rotation().to_matrix()

            self.integration_memos[key] = axial * bone_axes.inverted()
        return self.integration_memos[key]

    def rotate_by_world_vector(self, rot_v):
        local_rot_v = self.get_world_to_axial() * rot_v

        axis = local_rot_v.normalized()
        angle = local_rot_v.length

        rot = Quaternion(axis, angle)

        self.new_local_axial_rotation.rotate(rot)

    # Find new position and rotation based on force and torque
    def integrate(self, dt):
        # angular
        angular_accel = self.compute_angular_acceleration()
        self.new_angular_velocity += angular_accel * dt

        self.rotate_by_world_vector(self.new_angular_velocity * dt) 
        for child in self.get_child_bones():
            child.rotate_by_world_vector(-self.new_angular_velocity * dt)

        # linear
        self.new_linear_velocity += self.compute_linear_acceleration() * dt

        if self.bone.use_local_location():
            self.translation = self.get_world_to_axial() * self.compute_head_velocity() * dt
        else:
            self.translation = self.new_linear_velocity * dt

        self.integration_memos = [None, None, None, None, None]

    def apply(self):
        '''
        Apply updated values for rotation, angular velocity, and linear velocity
        to the wrapped bone
        '''
        self.bone.set_axial_rotation(self.new_local_axial_rotation)
        self.bone.set_angular_velocity(self.new_angular_velocity)
        self.bone.set_linear_velocity(self.new_linear_velocity)   

        if self.bone.get_parent() == None:
            self.bone.translate(self.translation)

class CPWrapper:
    '''
    Wrapper for the control point that allows for interpolation
    '''
    def __init__(self, parent, control_point):
        self.parent = parent
        self.control_point = control_point

        # self.interpolated_position = self.control_point.get_position()

    @property
    def interpolated_position(self):
        key = 'fwik_cp_ipos'
        if key not in self.control_point.object:
            self.control_point.object[key] = self.get_bone().get_tail_position()
        return Vector(self.control_point.object[key])
            
    @interpolated_position.setter
    def interpolated_position(self, value):
        key = 'fwik_cp_ipos'
        self.control_point.object[key] = value

    def get_bone(self):
        return self.parent.get_physbone( self.control_point.get_bone() )

    def advance_position(self):
        last_pos = self.interpolated_position
        next_pos = self.control_point.get_position() 

        max_iter = self.parent.iterations * 0.25
        iteration = self.parent.iteration

        if iteration >= max_iter:
            self.interpolated_position = next_pos
            return

        alpha = 1.0 / (max_iter - iteration)

        def linear_interpolation():
            self.interpolated_position = (1 - alpha) * last_pos + alpha * next_pos

        # Set origin of rotation to head of attached bone
        # origin = self.get_bone().get_head_position()

        # Set origin of rotation to center of mass
        origin = self.parent.compute_center_of_mass()

        def orbit():
            if (last_pos - origin).length < 0.001 or (next_pos - origin).length < 0.001:
                linear_interpolation()
                return

            # Rotate around origin
            last_r = (last_pos - origin).length
            next_r = (next_pos - origin).length
            new_r = (1 - alpha) * last_r + alpha * next_r

            dot = (last_pos - origin).dot(next_pos - origin) / (last_r * next_r)
            if dot > 1:
                dot = 1
            if dot < -1:
                dot = -1

            angle = math.acos(dot)

            i = (last_pos - origin).normalized()
            k = i.cross(next_pos - origin).normalized()
            j = k.cross(i).normalized()

            self.interpolated_position = origin + math.cos(angle * alpha) * new_r * i + math.sin(angle * alpha) * new_r * j

        # orbit()
        linear_interpolation()

    def get_position(self):
        return self.interpolated_position
        # return self.control_point.get_position()

class Simulator:
    '''
    Runs a FWIK simulation on the rig

    Keyword Arguments

    rig -- The rig to run the simulation on. Must have the following properties:
        - get_control_points() - returns a list of control points. See more below.
        - get_bones() - returns a list of bones. See more below.
        - damping - the amount of damping to apply in the simulation
    iterations -- number of iterations to run the simulation for
    time_step -- the time step each the iteration
    after_step -- a function to call after each step
    iteration -- (For debugging) run the simulation at certain iteration (useful for
            debugging control point interpolation

    In addition, the control points in rig.get_control_points() and the bones in
    rig,get_bones() should have the following properties

    ControlPoint:
        - get_position()
        - get_bone()
        - get_attachment_position()
        - get_spring_constant()

    Bone:
        - get_key() - returns some unique identifier that won't change during a
                    single run of the simulation
        - get_parent() - gets the parent of this bone. None if the bone is root

        - get_average_radius() - gets average radius of the bone

        - get_head_position() - gets the position of the head of the bone
        - get_tail_position() - gets the position of the tail of the bone
        - get_length() - gets the length of the bone
        - get_axial_rotation() - gets the rotation of the bone relative to its
                                parent (in quaternions)
        - set_axial_rotation() - sets the rotation of the bone relative to its
                                parent (in quaternions)
        - get_world_rotation() - gets the rotation of the bone in world space
        - get_mass() - gets the mass of the bone

        - get_angular_velocity() - gets the initial angular velocity of the bone 
                                (used for testing dynamics)
        - set_angular_velocity() - sets the initial angular velocity of the bone
                                for the next run of the simulation (used for 
                                testing dynamics)
        - get_linear_velocity() - gets the linear angular velocity of the bone's
                                COM (used for testing dynamics)
        - set_linear_velocity() - sets the initial angular velocity of the bone's
                                COM for the next run of the simulation (used for 
                                testing dynamics)

        

        * Bones are assumed to rotate at the head, and only the root bone is assumed
        to be translatable
    '''

    def __init__(self, rig, iterations=50, time_step=0.05, after_step=None, iteration=None):
        self.rig = rig
        self.iterations = iterations
        self.time_step = time_step
        self.after_step = after_step

        self.bone_map = dict()

        self.iteration = iteration

        self.root_bone = None

        self.joints = None

        self.center_of_mass = None

    def get_physbone(self, bone):
        if bone == None:
            return None
        else:
            key = bone.get_key()
            if key in self.bone_map:
                return self.bone_map[key]
            else:
                self.bone_map[key] = PhysBone(self, bone)
                return self.bone_map[key]

    def get_control_point(self, cp):
        return CPWrapper(self, cp)

    # def get_root_bone(self):
    #     if self.root_bone != None:
    #         return self.root_bone
    #     else:
    #         self.root_bone = self.get_physbone(self.rig.get_bones()[0])
    #         parent_bone = self.root_bone.get_parent_bone()
    #         while parent_bone != None:
    #             self.root_bone = parent_bone
    #             parent_bone = self.root_bone.get_parent_bone()
    #         return self.root_bone

    def compute_center_of_mass(self):
        if self.center_of_mass == None:
            bones = [self.get_physbone(bone) for bone in self.rig.get_bones()]

            com = Vector([0,0,0])
            total_mass = 0
            for bone in bones:
                com += bone.get_center_position() * bone.mass
                total_mass += bone.mass
            self.center_of_mass = com / total_mass
        return self.center_of_mass

    def get_joints(self):
        if self.joints == None:
            bones = [self.get_physbone(bone) for bone in self.rig.get_bones()]
            self.joints = []
            visited = set()
            for bone in bones:
                parent = bone.get_parent_bone()
                if parent != None and parent not in visited:
                    visited.add(parent)
                    neighbors = [parent] + parent.get_child_bones()
                    self.joints.append(neighbors)
        return self.joints

    def run(self):
        bones = [self.get_physbone(bone) for bone in self.rig.get_bones()]
        cps = [self.get_control_point(cp) for cp in self.rig.get_control_points()]

        dt = self.time_step

        def step(iteration):
            def torsion_forces(bone):
                tolerance = 0.1
                k = 10

                parent = bone.get_parent_bone()

                # Compute spring bounds
                min_rx = bone.bone.get_min_rot_x()
                max_rx = bone.bone.get_max_rot_x()
                range_rx = max_rx - min_rx

                min_ry = bone.bone.get_min_rot_y()
                max_ry = bone.bone.get_max_rot_y()
                range_ry = max_ry - min_ry

                min_rz = bone.bone.get_min_rot_z()
                max_rz = bone.bone.get_max_rot_z()
                range_rz = max_rz - min_rz

                # Compute rotations
                rot = bone.new_local_axial_rotation.to_exponential_map()

                def residue(angle):
                    angle = angle % (2 * math.pi)
                    if angle > math.pi:
                        angle -= 2 * math.pi
                    if angle < -math.pi:
                        angle += 2 * math.pi
                    return angle

                rot.x = residue(rot.x)
                rot.y = residue(rot.y)
                rot.z = residue(rot.z)

                ov_rot_x = rot.x - (max_rx - range_rx * tolerance)
                un_rot_x = rot.x - (min_rx + range_rx * tolerance)

                ov_rot_y = rot.y - (max_ry - range_ry * tolerance)
                un_rot_y = rot.y - (min_ry + range_ry * tolerance)

                ov_rot_z = rot.z - (max_rz - range_rz * tolerance)
                un_rot_z = rot.z - (min_rz + range_rz * tolerance)

                # Compute torque (linear w.r.t. angle)
                delta_x = 0
                delta_y = 0
                delta_z = 0

                if bone.bone.use_ik_limit_x():
                    if ov_rot_x > 0:
                        delta_x = ov_rot_x
                    elif un_rot_x < 0:
                        delta_x = un_rot_x
                
                if bone.bone.use_ik_limit_y():
                    if ov_rot_y > 0:
                        delta_y = ov_rot_y
                    elif un_rot_y < 0:
                        delta_y = un_rot_y
                
                if bone.bone.use_ik_limit_z():
                    if ov_rot_z > 0:
                        delta_z = ov_rot_z
                    elif un_rot_z < 0:
                        delta_z = un_rot_z
                
                delta = Vector([delta_x, delta_y, delta_z])

                #linear
                torque = -k * delta
                # capped linear
                # torque = -k * delta / max(delta.length, 1)

                # Apply torque
                world_torque = bone.get_world_to_axial().inverted() * torque
                bone.apply_torque(world_torque)
                if parent:
                    parent.apply_torque(-world_torque)
            
            def internal_forces():
                # max_diff = 0
                # for bone in bones:
                #     parent = bone.get_parent_bone()
                #     if parent == None:
                #         # M = bone.compute_force_acceleration_matrix(bone.get_head_position())
                #         # ha = bone.compute_head_acceleration()

                #         # diff = ha

                #         # Fr = M.inverted() * -ha
                #         # bone.apply_force(Fr, bone.get_head_position())
                #         pass
                #     else:
                #         M1 = parent.compute_force_acceleration_matrix(parent.get_tail_position())
                #         M2 = bone.compute_force_acceleration_matrix(bone.get_head_position())
                #         pta = parent.compute_tail_acceleration()
                #         cha = bone.compute_head_acceleration()

                #         diff = pta - cha

                #         Fr = (M1 + M2).inverted() * (cha - pta)
                #         parent.apply_force(Fr, parent.get_tail_position())
                #         bone.apply_force(-Fr, bone.get_head_position())

                #         max_diff = max(max_diff, abs(diff.length)))

                # return max_diff

                max_diff = 0
                for neighbors in self.get_joints():
                    accelerations = [neighbors[0].compute_tail_acceleration()] + [child.compute_head_acceleration() for child in neighbors[1:]]
                    diff = max((a - b).length for a, b in combinations( accelerations, 2 ))
                    max_diff = max(max_diff, diff)
                    
                    A = [None for _ in range(len(neighbors))]

                    A[0] = neighbors[0].compute_force_acceleration_matrix( neighbors[0].get_tail_position() ).to_4x4()
                    A[0] = (Matrix.Translation( neighbors[0].compute_tail_acceleration() ) * A[0]).inverted()

                    for i in range(1, len(neighbors)):
                        A[i] = neighbors[i].compute_force_acceleration_matrix( neighbors[i].get_head_position() ).to_4x4()
                        A[i] = (Matrix.Translation( neighbors[i].compute_head_acceleration() ) * A[i]).inverted()

                    v = sum(A[1:], A[0]).inverted() * Vector([0,0,0,1])
                    v = v / v.w
                    
                    neighbors[0].apply_force( (A[0] * v).xyz, neighbors[0].get_tail_position() )
                    for i in range(1, len(neighbors)):
                        neighbors[i].apply_force( (A[i] * v).xyz , neighbors[i].get_head_position() )

                return max_diff

            # print("Iteration {}".format(iteration))

            # Body forces and damping
            for bone in bones:
                bone.reset()
                # Linear damping
                bone.apply_force(-self.rig.damping * bone.compute_head_velocity(), bone.get_head_position())
                bone.apply_force(-self.rig.damping * bone.compute_tail_velocity(), bone.get_tail_position())
                # bone.apply_force(-self.rig.damping * bone.new_linear_velocity, bone.get_center_position())
                # Rotational damping
                # bone.apply_torque(-self.rig.damping * bone.new_angular_velocity)
                bone.new_angular_velocity *= 0.5

            # External Forces (control springs)
            for cp in cps:
                if cp.control_point.active:
                    cp.advance_position()
                    bone = cp.get_bone()

                    # Vector from control point to attachment point
                    r = cp.get_position() - cp.control_point.get_attachment_position()

                    # Hooke's law / Linear spring force with rest length = 0
                    # f = cp.control_point.get_spring_constant() * r

                    # Capped linear force spring (constant/linear hybrid)
                    f = cp.control_point.get_spring_constant() * r / max(r.length, 1)

                    bone.apply_force(f, cp.control_point.get_attachment_position())

            # Torsion forces
            for bone in bones:
                torsion_forces(bone)

            # Internal Forces (joint forces)
            iters = 0
            for _ in range(5 * len(bones) * len(bones)):
                iters += 1
                max_diff = internal_forces()
                if max_diff < 0.001:
                    break
            # print(iters, max_diff)

            # integrate
            for bone in bones:
                bone.integrate(dt) 
                bone.apply()

                self.center_of_mass = None

            # Use Blender's FK instead writing our own
            bpy.context.scene.update()

            if self.after_step != None:
                self.after_step()

        if self.iteration != None:
            step(self.iteration)
        else:
            for iteration in range(self.iterations):
                self.iteration = iteration
                step(iteration)
