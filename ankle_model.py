import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from data import Data

class FootDropAnkleModel:
  
  # constructor, put constants here
    def __init__(self):
        self.foot_mass = 1.0275  # kg
        self.centre_of_foot_mass_wrt_ankle = 0.1145 # cf -> m --> original 11.45  # cm
        self.interia_of_foot_around_ankle = 0.0197  # J -> (kg*m^2) 
        self.ta_moment_arm_wrt_ankle = 0.037 # 3.7[cm]  # d -> m
            
        self.time_activation = 0.01  # seconds
        self.time_deactivation = 0.04  # seconds

        self.viscosity_parameter = 0.82  # no units

        self.av_param = 1.33  # av-> no units, first force velocity paramter
        self.fv1 = 0.18  # fv1 -> no units, 2nd force velocity param
        self.fv2 = 0.023  # fv2 -> no units, 3rd force velocity parameter
        self.vmax = -0.9  # vmax -> m/s (shortening)
        self.max_isometric_force = 600  # N
        
        # defines range of muscle displacement where a force still remains perceptible
        self.shape_param = 0.56  # shape param of f_fl
        self.const_tendon_length = 0.223 # 22.3 [cm] lt -> m
        self.muscle_tendon_length_rest = 0.321 # m --> original 32.1 # lmt,0 -> cm
        self.elastic_torque_params = [2.1, -0.08, -7.97, 0.19, -1.79] # [a1, a2, a3, a4], Tela
        self.muscle_excitation_vector = Data.muscle_excitation_level_fig3
    #     self.ext1 = Data.linear_acc_shank_x
    #     self.ext2 = Data.linear_acc_shank_z
    #     self.ext3 = Data.abs_orientation_shank
    #     self.ext4 = Data.abs_velocity_rotation_shank
        self.external_vector = [Data.linear_acc_shank_x, Data.linear_acc_shank_z, 
                                Data.abs_orientation_shank, Data.abs_velocity_rotation_shank]


    def normalize_times(self):
    #     self.ext1[:,0] /= 0.54
    #     self.ext2[:,0] /= 0.54
    #     self.ext3[:,0] /= 0.54
    #     self.ext4[:,0] /= 0.54
        for state in self.external_vector:
            state[:,0] /= 0.54
        self.external_vector[2][:,1] /= 5
        self.external_vector[3][:,1] /= 5
        self.muscle_excitation_vector[:,0] /= 0.35
        self.muscle_excitation_vector[:,] += 1


    def get_current_muscle_excitation(self, time):
        for i in range(len(self.muscle_excitation_vector)):
            if self.muscle_excitation_vector[i, 0] > time:
                return self.muscle_excitation_vector[i, 1]


    # def get_current_ext_vector(self, x_ext, t):
    #     """
    #     :param t: time
    #     :param x_ext: 2D array containing time poins and data points for external state vecctors
    #     :return data value of external state vector closest to passed in time point: 
    #     """
    #     # iterate through all data points in corresponding data array
    #     # swing phase starts for external state vector at 0.55s, so time is adjusted to this point 
    #     # in simulation
    #     norm_x_ext = x_ext
    #     norm_x_ext[:,0] /= 0.54
    #     for i in range(len(norm_x_ext)-1):
    #         if norm_x_ext[i, 0] > t:
    #             return norm_x_ext[i-1, 1]

    def get_current_ext_vector(self, time):
        current_values = []
        # for state in range(len(self.external_vector)):
        for state in self.external_vector:
            for i in range(len(state)-1):
                if state[i,0] > time:
                    # if state.all() == self.external_vector[2].all():
                    #     state[i-1,1] *= -1
                    current_values.append(state[i-1,1])
            # for i in range(len(self.external_vector[state])-1):
            #     if self.external_vector[state][i,0] > time:
            #         if i == 2:
            #             self.external_vector[state][i-1, 1] *= -1
            #         current_values.append(self.external_vector[state][i-1,1])
                    break
        return current_values


    def compute_gravity_torque(self, x2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :return: gravity torque of the foot around the ankle
        """
        g = 9.81  # acceleration of gravity m/s^2
        return -(self.foot_mass * self.centre_of_foot_mass_wrt_ankle * math.cos(x2) * g)


    def compute_acceleration_torque(self, x2, x_ext1, x_ext2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :param x_ext1: linear acceleration of foot in horizontal direction
        :param x_ext2: linear acceleration of foot in vertical direction
        :return: torque induced by the movement of the ankle
        """
        return self.foot_mass * self.centre_of_foot_mass_wrt_ankle * (x_ext1 * math.sin(x2) - 
               x_ext2 * math.cos(x2))


    def compute_force_velocity(self, x3, x_ext4):
        """
        :param t: time in seconds
        :param x3: abs rotational velocity
        :param x_ext4: absolute velocity rotation of the shank
        :return: force velocity
        """
        
        velocity_ce = self.ta_moment_arm_wrt_ankle * (x_ext4 - x3)

        # velocity_ce < 0 -> contraction
        if velocity_ce < 0:
            return (1 - (velocity_ce/self.vmax))/(1 + (velocity_ce/(self.vmax*self.fv1)))
        
        # velocity_ce >= 0 -> extension or isometric
        else:
            return (1 + (self.av_param*(velocity_ce/self.fv2)))/(1 + (velocity_ce/self.fv2))


    def compute_ta_muscular_force(self, x1, x2, x3, x_ext3, x_ext4):
        """
        :param x1: dyanimc activation level of foot
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :param x3: abs rotational velocity
        :param x_ext3: absolute rotation orientation of the shank
        :param x_ext4: absolute velocity rotation of the shank
        :return: tibialus anterior muscular force
        """
        # non-linear relationship linking generated force to length of muscle and therefore 
        # to ankle joint angle
        length_muscle_tendon = self.muscle_tendon_length_rest - self.ta_moment_arm_wrt_ankle * (x_ext3 - x2)
        length_ce =  length_muscle_tendon - self.const_tendon_length

        length_ce_opt = 0.321 # m
        f_fl_eqn = (length_ce - length_ce_opt) / (self.shape_param*length_ce_opt)
        f_fl = math.exp(-math.pow(f_fl_eqn,2))
        f_fv = self.compute_force_velocity(x3, x_ext4)

        return x1 * self.max_isometric_force * f_fl * f_fv 


    def compute_passive_elastic_torque(self, x2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :return: passive elastic torque around the ankle
        :note: it is assumed that the knee angle variation is not significant to the ankle angular 
        range 
        """
        a1 = self.elastic_torque_params[0]
        a2 = self.elastic_torque_params[1]
        a3 = self.elastic_torque_params[2]
        a4 = self.elastic_torque_params[3]
        a5 = self.elastic_torque_params[4]

        return math.exp(a1+a2*x2) - math.exp(a3+a4*x2) + a5


    def compute_state_vector_derivative(self, t, x):
        """
        :param t: time
        :param x: state vector [dyanimc activation level of foot, abs orientation of foot, abs 
        rotational velocity of foot]
        :return: derivative of state vector
        """
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]

        # external state vector
        # x_ext1 = self.get_current_ext_vector(Data.linear_acc_shank_x, t)
        # x_ext2 = self.get_current_ext_vector(Data.linear_acc_shank_z, t)
        # x_ext3 = self.get_current_ext_vector(Data.abs_orientation_shank, t) 
        # x_ext4 = self.get_current_ext_vector(Data.abs_velocity_rotation_shank, t) 
        current_external_vector = self.get_current_ext_vector(t)
        x_ext1 = current_external_vector[0]
        x_ext2 = current_external_vector[1]
        x_ext3 = current_external_vector[2] 
        x_ext4 = current_external_vector[3]
        # get parameters needed to calculate derivatives
        Fm = self.compute_ta_muscular_force(x1, x2, x3, x_ext3, x_ext4)
        T_eta = self.compute_passive_elastic_torque(x2)
        T_grav = self.compute_gravity_torque(x2)
        T_acc = self.compute_acceleration_torque(x2, x_ext1, x_ext2)
        muscle_excitation = self.get_current_muscle_excitation(t)
        
        # calc derivative of state vector
        dt_x1 = (muscle_excitation - x1) * ((muscle_excitation / self.time_activation) - ((1 - muscle_excitation) / self.time_deactivation)) # (unitless activation)/seconds
        dt_x2 = float(x3) # degrees/s
        dt_x3 = (1/self.interia_of_foot_around_ankle)*(Fm*self.ta_moment_arm_wrt_ankle + T_grav + T_acc + T_eta + self.viscosity_parameter*(x_ext4-x3)) # degrees/s^2
        
        # return derivative state vector
        return np.array([dt_x1, dt_x2, dt_x3])


    def simulate(self, sim_duration):
        """
        :param sim_duration: duration of the simulation in seconds
        :return: solution to state derivative IVP
        """
        # normalize time
        x_init = [0.25, -9, 1]
        sol = scipy.integrate.solve_ivp(self.compute_state_vector_derivative, [1, 2], x_init, method='RK45', max_step=0.01)
        return sol