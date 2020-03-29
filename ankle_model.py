import numpy as np
import math

class FootDropAnkleModel:
  
  # constructor, put constants here
  def __init__(self):
    self.foot_mass = 1.0275  # kg
    self.centre_of_foot_mass_wrt_ankle =11.45  # cm
    self.interia_of_foot_around_ankle = 0.0197  # J (kg*m^2) 
    self.ta_moment_arm_wrt_ankle = 3.7  # cm
        
    self.time_activation = 0.01  # seconds
    self.time_deactivation = 0.04  # seconds

    self.viscosity_parameter = 0.82  # no units

    self.av_param = 1.33  # no units, first force velocity paramter
    self.fv1 = 0.18  # no units, 2nd force velocity param
    self.fv2 = 0.023  # no units, 3rd force velocity parameter
    self.max_contraction_speed = -0.9  # m/s (shortening)
    self.max_isometric_force = 600  # N
    # defines range of muslce displacement where a force still remains perceptible
    self.shape_param = 0.56  # shape param of f_fl
    self.const_tendon_length = 22.3 # cm
    self.muscle_tendon_length_rest = 32.1 # cm
    self.elastic_torque_params = [2.1, -0.08, -7.97, 0.19, -1.79] # [a1, a2, a3, a4], Tela

    self.muscle_excitation = [0, 1] # TODO: figure out how to calculate this or initialize!!!!!!!!!!!!!!!! AKA: u1 in the paper

  # state vector x = [x1, x2, x3] = [fact, alpha_f, rot_vel_f]
  # x1 --> Dynamic activation level of foot, values 0 to 1
  # derivative of state vector: long af
	
  # external state vector function defs: 
  # linear acceleration in horizontal and vertical directions
  # absolute velocity rotation of the shank
  # ^^ THIS DATA NEEDS TO BE DONE IN REAL TIME --> We will take it from the graphs and use their experiemental data
  
  
  #eqn 4
    def compute_gravity_torque(self, alpha):
    """
    :param alpha: absolute orientation of foot wrt horizontal axis 
    :return: gravity torque of the foot around the ankle
    """
        g = 9.81  # acceleration of gravity
        return −(self.foot_mass * self.centre_of_foot_mass_wrt_ankle * math.cos(alpha) * g)


    # eqn 5
    def compute_acceleration_torque(self):
    """
    :return: torque induced by the movement of the ankle
    """

    # eqn 6
    def compute_elastic_torque(self):
    
    
  # eqn 7
    def compute_ta_muscular_force(self):
        # non-linear relationship linking generated force to length of muscle and therefore to ankle joint angle
        length_tendon = self.const_tendon_length
        mt_length_rest = self.muscle_tendon_length_rest
        length_muscle_tendon = mt_length_rest + self.ta_moment_arm_wrt_ankle(x3_ext - x2)
        # length_ce =  # lmt - lt
        f_fl = math.exp()
    
    
  # state variable 
    def dynamics(self, x):
    """
    :param x: state vector 
    :return: derivative of state vector
    """
    
    # TODO: get parameters needed to calculate derivatives
        Fm = self.compute_ta_muscular_force() # TODO fill parameters
        T_grav = self.compute_gravity_torque(x[2])
        T_acc = self.compute_acceleration_torque() # TODO fill parameters

        dt_x1 = (self.muscle_exitation - x[1]) * ((self.muscle_exitation / self.time_activation) - ((1 - self.muscle_exitation) / self.time_deactivation))
        dt_x2 = x[3]
        dt_x3 = 1/self.interia_of_foot_around_ankle(self) + T_grav + T_acc
  
