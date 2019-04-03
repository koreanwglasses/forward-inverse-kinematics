import bpy
from mathutils import *

class PhysBone:
    '''
    A wrapper for the bone to perform calculations and modify the actual bone
    '''
    def __init__(self, parent, bone):
        self.parent = parent
        self.bone = bone

        # Physics parameters
        self.new_axial_rotation = self.bone.get_axial_rotation()
        self.new_angular_velocity = self.bone.get_angular_velocity()
        self.new_linear_velocity = self.bone.get_linear_velocity()

        self.net_force = Vector((0,0,0))
        self.net_torque = Vector((0,0,0))

        # Compute principle moments asusming bone is a cylinder
        m = self.bone.get_mass()
        r = self.bone.get_average_radius()
        h = self.bone.get_length()
        principle_moment_x = (1.0 / 12.0) * m * (3 * r**2 + h**2)
        principle_moment_y = (1.0 / 2.0) * m * r**2
        self.principle_moment = Matrix([[principle_moment_x, 0, 0 ], 
                                    [0, principle_moment_y, 0],
                                    [0, 0, principle_moment_x]])

        self.translation = Vector([0,0,0])

    def get_parent_bone(self):
        return self.parent.get_physbone(self.bone.get_parent())

    def get_child_bone(self):
        return self.parent.get_physbone(self.bone.get_child())

    def get_head_position(self):
        return self.bone.get_head_position()

    def get_tail_position(self):
        return self.bone.get_tail_position()

    def get_center_position(self):
        return 0.5 * (self.bone.get_head_position() + self.bone.get_tail_position())
    
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
        Icm = self.compute_moment()
        return Icm.inverted() * (self.net_torque - self.new_angular_velocity.cross(Icm * self.new_angular_velocity))

    def compute_linear_acceleration(self):
        return self.net_force / self.bone.get_mass()

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

        return Matrix.Scale(1.0 / self.bone.get_mass(), 3) + rx.transposed() * self.compute_moment().inverted() * rx

    def compute_moment(self):
        world_rotation = self.get_world_rotation().to_matrix()
        return world_rotation * self.principle_moment * world_rotation.inverted()

    def reset(self):
        self.net_force = Vector((0,0,0))
        self.net_torque = Vector((0,0,0))

    def apply_force(self, force, position):
        self.net_force += force

        r = position - self.get_center_position()
        torque = r.cross(force)

        self.net_torque += torque

    def apply_torque(self, torque):
        self.net_torque += torque

    def apply_axial_torque(self, axial_torque):
        self.net_torque += self.get_world_to_axial().inverted() * axial_torque

    def get_world_to_axial(self):
        bone_axes = self.bone.get_bone_axes_to_world()
        axial = self.bone.get_axial_rotation().to_matrix()

        return axial * bone_axes.inverted()

    # Find new position and rotation based on force and torque
    def integrate(self, dt):
        # angular
        angular_accel = self.compute_angular_acceleration()
        self.new_angular_velocity += angular_accel * dt

        axial_angular_velocity = self.get_world_to_axial() * self.new_angular_velocity

        axis = axial_angular_velocity.normalized()
        angle = axial_angular_velocity.length

        delta_rot = Quaternion(axis, angle * dt)

        self.new_axial_rotation.rotate(delta_rot)

        # linear
        self.new_linear_velocity += self.compute_linear_acceleration() * dt
        self.translation = self.compute_head_velocity() * dt
        # self.translation += self.new_linear_velocity * dt

    def apply(self):
        '''
        Apply updated values for rotation, angular velocity, and linear velocity
        to the wrapped bone
        '''
        self.bone.set_axial_rotation(self.new_axial_rotation)
        self.bone.set_angular_velocity(self.new_angular_velocity)
        self.bone.set_linear_velocity(self.new_linear_velocity)   

        if self.bone.get_parent() == None:
            self.bone.translate(self.translation)

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

    def __init__(self, rig, iterations=50, time_step=0.05, after_step=None):
        self.rig = rig
        self.iterations = iterations
        self.time_step = time_step
        self.after_step = after_step

        self.bone_map = dict()

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

    def run(self):
        bones = [self.get_physbone(bone) for bone in self.rig.get_bones()]
        cps = self.rig.get_control_points()

        dt = self.time_step

        def step():
            def torsion_forces(bone):
                tolerance = .1
                k = 20

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
                rot = bone.new_axial_rotation.to_exponential_map()

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
                bone.apply_axial_torque(torque)
                if parent:
                    parent.apply_axial_torque(-torque)


            def internal_forces():
                max_diff = 0
                for bone in bones:
                    parent = bone.get_parent_bone()
                    if parent == None:
                        # M = bone.compute_force_acceleration_matrix(bone.get_head_position())
                        # ha = bone.compute_head_acceleration()

                        # diff = ha

                        # Fr = M.inverted() * -ha
                        # bone.apply_force(Fr, bone.get_head_position())
                        pass
                    else:
                        M1 = parent.compute_force_acceleration_matrix(parent.get_tail_position())
                        M2 = bone.compute_force_acceleration_matrix(bone.get_head_position())
                        pta = parent.compute_tail_acceleration()
                        cha = bone.compute_head_acceleration()

                        diff = pta - cha

                        Fr = (M1 + M2).inverted() * (cha - pta)
                        parent.apply_force(Fr, parent.get_tail_position())
                        bone.apply_force(-Fr, bone.get_head_position())

                        max_diff = max(max_diff, abs(diff.length))
                return max_diff

            # print("~~~~~")

            # Body forces and damping
            for bone in bones:
                bone.reset()
                # Linear damping
                bone.apply_force(-self.rig.damping * bone.compute_head_velocity(), bone.get_head_position())
                bone.apply_force(-self.rig.damping * bone.compute_tail_velocity(), bone.get_tail_position())
                # bone.apply_force(-self.rig.damping * bone.new_linear_velocity, bone.get_center_position())
                # Rotational damping
                bone.apply_torque(-self.rig.damping * bone.new_angular_velocity)

                # print(bone.new_axial_rotation)

            # External Forces (control springs)
            for cp in cps:
                if cp.active:
                    bone = self.get_physbone(cp.get_bone())

                    # Vector from control point to attachment point
                    r = cp.get_position() - cp.get_attachment_position()

                    # Hooke's law / Linear spring force with rest length = 0
                    # f = cp.get_spring_constant() * r

                    # Capped linear force spring (constant/linear hybrid)
                    f = cp.get_spring_constant() * r / max(r.length, 1)

                    bone.apply_force(f, cp.get_attachment_position())

            # Torsion forces
            for bone in bones:
                torsion_forces(bone)

            # Internal Forces (joint forces)
            # iters = 0
            for _ in range(5 * len(bones) * len(bones)):
                # iters += 1
                max_diff = internal_forces()
                if max_diff < 0.001:
                    break
            # print(iters, max_diff)

            # integrate
            for bone in bones:
                bone.integrate(dt) 
                bone.apply()

            # Use Blender's FK instead writing our own
            bpy.context.scene.update()

            if self.after_step != None:
                self.after_step()

        for _ in range(self.iterations):
            step()