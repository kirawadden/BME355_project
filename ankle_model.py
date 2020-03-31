import numpy as np
import math
import scipy.integrate
from data import Data

class FootDropAnkleModel:
  
  # constructor, put constants here
    def __init__(self):
        self.foot_mass = 1.0275  # kg
        self.centre_of_foot_mass_wrt_ankle =11.45  # cm
        self.interia_of_foot_around_ankle = 0.0197  # J -> (kg*m^2) 
        self.ta_moment_arm_wrt_ankle = 3.7  # d -> cm
            
        self.time_activation = 0.01  # seconds
        self.time_deactivation = 0.04  # seconds

        self.viscosity_parameter = 0.82  # no units

        self.av_param = 1.33  # av-> no units, first force velocity paramter
        self.fv1 = 0.18  # fv1 -> no units, 2nd force velocity param
        self.fv2 = 0.023  # fv2 -> no units, 3rd force velocity parameter
        self.vmax = -0.9  # vmax -> m/s (shortening)
        self.max_isometric_force = 600  # N
        # defines range of muslce displacement where a force still remains perceptible
        self.shape_param = 0.56  # shape param of f_fl
        self.const_tendon_length = 22.3 # lt -> cm
        self.muscle_tendon_length_rest = 32.1 # lmt,0 -> cm
        self.elastic_torque_params = [2.1, -0.08, -7.97, 0.19, -1.79] # [a1, a2, a3, a4], Tela

        # TODO: MAKE THIS VECTOR NOT CONST
        self.muscle_excitation = 0.5 #[0, 1] # TODO: figure out how to calculate this or initialize!!!!!!!!!!!!!!!! AKA: u1 in the paper

  # state vector x = [x1, x2, x3] = [fact, alpha_f, rot_vel_f]
  # x1 --> Dynamic activation level of foot, values 0 to 1
  # derivative of state vector: long af
	
  # external state vector function defs: 
  # linear acceleration in horizontal and vertical directions
  # absolute velocity rotation of the shank
  # ^^ THIS DATA NEEDS TO BE DONE IN REAL TIME --> We will take it from the graphs and use their experiemental data
  
    def get_current_ext_vector(self, x_ext, t):
        """
        :param t: time
        :param x_ext: 2D array containing time poins and data points for external state vecctors
        :return data value of external state vector closest to passed in time point: 
        """
        # TODO: return value closer to time vs. just at i
        # TODO: ensure that this works properly
        # iterate through all data points in corresponding data array
        for i in range(len(x_ext)-1):
            # and x_ext[i+1][0] <= t
            if x_ext[i][0] > t:
                return x_ext[i-1][1]
    
  
  #eqn 4
    def compute_gravity_torque(self, alpha):
        """
        :param alpha: absolute orientation of foot wrt horizontal axis 
        :return: gravity torque of the foot around the ankle
        """
        g = 9.81  # acceleration of gravity
        return -(self.foot_mass * self.centre_of_foot_mass_wrt_ankle * math.cos(alpha) * g)


  # eqn 5
    def compute_acceleration_torque(self, x_ext1, x_ext2, alpha):
        """
        :return: torque induced by the movement of the ankle
        """
        return self.foot_mass*self.centre_of_foot_mass_wrt_ankle*(x_ext1*math.sin(alpha) - x_ext2*math.cos(alpha))


    def compute_force_velocity(self, t, x):
        """
        :param t: time
        :param x: state vector
        :return: force velocity
        """
        x_ext4 = self.get_current_ext_vector(Data.abs_velocity_rotation_shank, t)
        velocity_ce = self.ta_moment_arm_wrt_ankle * (x_ext4 - x[1])

        # velocity_ce < 0 -> contraction
        if velocity_ce < 0:
            return (1 - (velocity_ce/self.vmax))/(1 + (velocity_ce/(self.vmax*self.fv1)))
        
        # velocity_ce >= 0 -> extension or isometric
        else:
            return (1 + (self.av_param*(velocity_ce/self.fv2)))/(1 + (velocity_ce/self.fv2))


    # eqn 7
    def compute_ta_muscular_force(self, t, x):
        """
        :param t: time
        :param x: state vector
        :return: tibialus anterior muscular force
        """
        # non-linear relationship linking generated force to length of muscle and therefore to ankle joint angle
        x_ext3 = self.get_current_ext_vector(Data.abs_orientation_shank, t)
        length_muscle_tendon = self.muscle_tendon_length_rest - self.ta_moment_arm_wrt_ankle * (x_ext3 - x[1])
        length_ce =  length_muscle_tendon - self.const_tendon_length
        # # TODO: find wtf the optimal length muscle tendon is
        length_ce_opt = 32.1
        f_fl_eqn = -(length_ce - length_ce_opt) / (self.shape_param*length_ce_opt)
        f_fl = math.exp(math.pow(f_fl_eqn,2))
        f_fv = self.compute_force_velocity(t, x)
        x_ext3 = self.get_current_ext_vector(Data.abs_orientation_shank, t)
        x_ext4 = self.get_current_ext_vector(Data.abs_velocity_rotation_shank, t)
        return x[0] * self.max_isometric_force * f_fl * f_fv * (x_ext4 - x[1])


    # eqn 6
    def compute_passive_elastic_torque(self, alpha):
        """
        :param alpha: x2 in state vector --> abs orientation of foot
        :return: passive elastic torque around the ankle
        :note: it is assumed that the knee angle variation is not significant to the ankle angular range 
        """
        a1 = self.elastic_torque_params[0]
        a2 = self.elastic_torque_params[1]
        a3 = self.elastic_torque_params[2]
        a4 = self.elastic_torque_params[3]
        a5 = self.elastic_torque_params[4]

        return math.exp(a1+a2*alpha) - math.exp(a3+a4*alpha) + a5

    # state variable eqns 1, 2, 3
    def compute_state_vector_derivative(self, t, x):
        """
        :param x: state vector (dyanimc activation level of foot, abs orientation of foot, abs rotational velocity of foot)
        :param t: time
        :return: derivative of state vector
        """
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]

        # external state vector
        x_ext1 = self.get_current_ext_vector(Data.linear_acc_shank_x, t)
        x_ext2 = self.get_current_ext_vector(Data.linear_acc_shank_z, t)
        x_ext3 = self.get_current_ext_vector(Data.abs_orientation_shank, t)
        x_ext4 = self.get_current_ext_vector(Data.abs_velocity_rotation_shank, t)

        # TODO: get parameters needed to calculate derivatives
        Fm = self.compute_ta_muscular_force(t, x)
        T_eta = self.compute_passive_elastic_torque(x2)
        T_grav = self.compute_gravity_torque(x2)
        T_acc = self.compute_acceleration_torque(x_ext1, x_ext2, x2)

        # calc derivative of state vector
        dt_x1 = (self.muscle_excitation - x1) * ((self.muscle_excitation / self.time_activation) - ((1 - self.muscle_excitation) / self.time_deactivation))
        dt_x2 = x3
        dt_x3 = 1/self.interia_of_foot_around_ankle*(Fm*self.ta_moment_arm_wrt_ankle + T_grav + T_acc + T_eta + self.viscosity_parameter*(x_ext4-x3))

        # return derivative state vector
        return [dt_x1, dt_x2, dt_x3]

    def simulate(self, sim_duration):
        """
        :param sim_duration: duration of the simulation in seconds
        """
        x = [0, 0, 0]  # TODO: maybe pass these into simulation

        return scipy.integrate.solve_ivp(self.compute_state_vector_derivative, (1, sim_duration), x, method='RK45', max_step=0.01)
        
