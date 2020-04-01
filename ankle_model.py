import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate
import random
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
        self.const_tendon_length = 0.223 # 22.3 [cm]lt -> m
        self.muscle_tendon_length_rest = 0.321 # m --> original 32.1 # lmt,0 -> cm
        self.elastic_torque_params = [2.1, -0.08, -7.97, 0.19, -1.79] # [a1, a2, a3, a4], Tela
  
    def get_current_muscle_excitation(self, t):
        for i in range(len(Data.muscle_excitation_level_fig3)):
            if Data.muscle_excitation_level_fig3[i][0] > t:
                return Data.muscle_excitation_level_fig3[i][1]

    def get_current_ext_vector(self, x_ext, t):
        """
        :param t: time
        :param x_ext: 2D array containing time poins and data points for external state vecctors
        :return data value of external state vector closest to passed in time point: 
        """
        # TODO: return value closer to time vs. just at i
        # TODO: ensure that this works properly
        # iterate through all data points in corresponding data array
        t += 0.55 ### added 0.55 to time as we believe this is where the swing phase starts
                  ### there's a pretty drastic magnitude different if we reduce to t += 0.5
        for i in range(len(x_ext)-1):
            if x_ext[i][0] > t:
                return x_ext[i-1][1]
    
  
  #eqn 4
    def compute_gravity_torque(self, x2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :return: gravity torque of the foot around the ankle
        """
        g = 9.81  # acceleration of gravity m/s^2
        return -(self.foot_mass * self.centre_of_foot_mass_wrt_ankle * math.cos(x2) * g)


  # eqn 5
    def compute_acceleration_torque(self, x2, x_ext1, x_ext2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :param x_ext1: linear acceleration of foot in horizontal direction
        :param x_ext2: linear acceleration of foot in vertical direction
        :return: torque induced by the movement of the ankle
        """
        return self.foot_mass * self.centre_of_foot_mass_wrt_ankle * (x_ext1 * math.sin(x2) - x_ext2 * math.cos(x2))


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


    # eqn 7
    def compute_ta_muscular_force(self, x1, x2, x3, x_ext3, x_ext4):
        """
        :param x1: dyanimc activation level of foot
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :param x3: abs rotational velocity
        :param x_ext3: absolute rotation orientation of the shank
        :param x_ext4: absolute velocity rotation of the shank
        :return: tibialus anterior muscular force
        """
        # non-linear relationship linking generated force to length of muscle and therefore to ankle joint angle
        length_muscle_tendon = self.muscle_tendon_length_rest - self.ta_moment_arm_wrt_ankle * (x_ext3 - x2)
        length_ce =  length_muscle_tendon - self.const_tendon_length

        length_ce_opt = 0.321 # 32.1 [cm] -> m
        f_fl_eqn = (length_ce - length_ce_opt) / (self.shape_param*length_ce_opt)
        f_fl = math.exp(-math.pow(f_fl_eqn,2))
        f_fv = self.compute_force_velocity(x3, x_ext4)

        return x1 * self.max_isometric_force * f_fl * f_fv 


    # eqn 6
    def compute_passive_elastic_torque(self, x2):
        """
        :param x2: alpha, absolute orientation of foot wrt horizontal axis 
        :return: passive elastic torque around the ankle
        :note: it is assumed that the knee angle variation is not significant to the ankle angular range 
        """
        a1 = self.elastic_torque_params[0]
        a2 = self.elastic_torque_params[1]
        a3 = self.elastic_torque_params[2]
        a4 = self.elastic_torque_params[3]
        a5 = self.elastic_torque_params[4]

        return math.exp(a1+a2*x2) - math.exp(a3+a4*x2) + a5

    # state variable eqns 1, 2, 3
    def compute_state_vector_derivative(self, t, x):
        """
        :param t: time
        :param x: state vector [dyanimc activation level of foot, abs orientation of foot, abs rotational velocity of foot]
        :return: derivative of state vector
        """
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]

        # external state vector
        x_ext1 = self.get_current_ext_vector(Data.linear_acc_shank_x, t)
        x_ext2 = self.get_current_ext_vector(Data.linear_acc_shank_z, t)
        x_ext3 = self.get_current_ext_vector(Data.abs_orientation_shank, t) #* (math.pi/180)
        x_ext4 = self.get_current_ext_vector(Data.abs_velocity_rotation_shank, t) # * (180/math.pi) ### might need unit conversion but also makes the graphs look like they resonate...

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
        x_init = [0, -15, 0]  
        sol = scipy.integrate.solve_ivp(self.compute_state_vector_derivative, [0.01, sim_duration], x_init, method='RK45', max_step = 0.01)
        return sol

    """
    Cannot start the simulation at 0. This is because the u input gets caught on the first value if we start at 0.
    The u vector data is only for the swing phase (.36 seconds). This is from the main paper.
    The main paper abs. orientation of shank ranges -15 to about +15 Figure 2
    The data paper abs. orientation of shank ranges -15 to about +50 Figure 6 a   
    Note, looking at data paper, abs. orientation of shank increases with stimulation. Normal walking has greatest range. 
    Therefore, we expect differences, seeing that the input data and euqations have different outputs.  *** angle orient. different expected values
    Current biggest concern: There is a big spike at the beginning of msucle activation (fact) - fig 3 graph c).
    We currently are not getting this big spike when we run the simulation. We think that the activation trend for our sim
    should following the exciation input vector.
    We think that it will be in our best interest if everyone can take a look at the code
    Thing we can look at:
    - Units (specifically radians vs degrees)
    - Double check all equations and constants
    - Possibly look at time step for equation solver
    - Check input data and verify that it matches the graphs that we extracted it fro
    - Verify that the methods being used to select the 'current' inputs (u and x_ext values) is correct
    Notes:
    - code looks difficult to read, but actually isn't too bad
    - Graphing happens in main.py, class methods are in ankle_model and data is in data
    - We are using muscle_excitation_level_fig3
    - Run the project by going to the file that you are storing this project in (in your terminal). Then type 'python3 main.py' into your terminal/command line.
    - We are running the sim for 0.35s as that is what they do in the paper, and the time length of the input vector (u).
        - We need to ensure that the data we are using from the external vector graphs corresponds to the swing phase (I think it start around 0.55s for those graphs but please feel free to double check this)
    
    Damage control:
    - Why didn't this work: using data from multiple papers, unclear what their ICs were, some constants undefined
    - How are we validating? Need to find other data to validate curve shape against? Maybe what they got is actually wack and what we are getting actually make somes sense
    - Double check figures/diagrams to ensure we are understanding all angles correctly
    - We can also talk about the range fo the graph - ex. 0 <= activation <= 1, which means that the values we are receving technically make sense, but the shape of the curve is off
        - Need to prove our model is somewhat realisitc :) 
    """