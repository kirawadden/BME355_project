import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from data import Data
from ankle_model import FootDropAnkleModel



def run_simulation():
    """
    Runs simulation for 0.35 s (swing phase). Plots solution (state vector).
    """
    foot_drop = FootDropAnkleModel()
    soln = foot_drop.simulate(0.35)
    plt.plot(soln.t, soln.y[0], 'r')
    plt.title('Muscle Activation vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Activation Level')
    plt.show()
    plt.plot(soln.t, soln.y[1], 'b')
    plt.title('Absolute Orientation of Foot vs. Time')
    plt.ylabel('Orientation (deg)')
    plt.xlabel('Time (s)')
    plt.show()
    plt.plot(soln.t, soln.y[2], 'g')
    plt.title('Absolute Rotational Velocity of Foot vs. Time')
    plt.ylabel('Velocity (deg/s)')
    plt.xlabel('Time (s)')
    plt.show()
    print(soln.message)
    print(soln.status)

def graph_data():
    """
    Graphs external state vector quantities and muscle excitation level. Can be used as a sanity
    check to visually verify that data extracted from web plot digitizer matches the original 
    graphs.
    """
    x_ext1 = np.array(Data.linear_acc_shank_x)
    x_ext2 = np.array(Data.linear_acc_shank_z)
    x_ext3 = np.array(Data.abs_orientation_shank) 
    x_ext4 = np.array(Data.abs_velocity_rotation_shank)
    u1 = np.array(Data.muscle_excitation_level_fig3)
    plt.plot(u1[:, 0], u1[:, 1])
    plt.title('Muscle Excitation vs. Time')
    plt.ylabel('Muslce Excitation Level')
    plt.xlabel('Time (s)')
    plt.show()
    plt.plot(x_ext1[:, 0], x_ext1[:, 1], label='x-direction')
    plt.plot(x_ext2[:, 0], x_ext2[:, 1], label='z-direction')
    plt.title('Linear Acceleration of the Shank vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.legend()
    plt.show()
    plt.plot(x_ext4[:, 0], x_ext4[:, 1])
    plt.title('Absolute Rotational Velocity of Shank vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Rotational Velocity (rad/s)')
    plt.show()
    plt.plot(x_ext3[:, 0], x_ext3[:, 1])
    plt.title('Absolute Orientation of Shank vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Orientation (deg)')
    plt.show()

def main():
    """
    Program starting point
    """
    run_simulation()
    graph_data()


if __name__ == '__main__':
    main()
