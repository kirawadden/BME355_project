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
    foot_drop.normalize_times()
    soln = foot_drop.simulate(0.35)
    plt.plot((soln.t-1)*0.35, soln.y[0]/2, 'r')
    plt.title('Muscle Activation vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Activation Level')
    plt.show()
    plt.plot((soln.t-1)*0.35, soln.y[1], 'b')
    plt.title('Orientation of Foot vs. Time')
    plt.ylabel('Orientation (deg)')
    plt.xlabel('Time (s)')
    plt.show()
    plt.plot((soln.t-1)*0.35, soln.y[2], 'g')
    plt.title('Rotational Velocity of Foot vs. Time')
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
    # external state vector 1 -> linear acceleration of shank x-direction
    x_ext1 = Data.linear_acc_shank_x
    # external state vector 2 -> linear acceleration of shank z-direction
    x_ext2 = Data.linear_acc_shank_z
    # external state vector 3 -> orientation of shank 
    x_ext3 = Data.abs_orientation_shank
    # external state vector 4 -> rotational velocity of shank
    x_ext4 = Data.abs_velocity_rotation_shank
    # input vector (u)
    u1 = Data.muscle_excitation_level_fig3

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
    plt.title('Rotational Velocity of Shank vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Rotational Velocity (rad/s)')
    plt.legend()
    plt.show()
    plt.plot(x_ext4[:, 0], x_ext4[:, 1])
    # partial swing phase
    plt.plot(x_ext4[6:15,0], x_ext4[6:15,1],'r', label='Swing Phase')
    plt.title('Rotational Velocity of Shank vs. Time (Partial Swing Phase)')
    plt.xlabel('Time (s)')
    plt.ylabel('Rotational Velocity (rad/s)')
    plt.legend()
    plt.show()
    plt.plot(x_ext4[:, 0], x_ext4[:, 1])
    # full swing phase
    plt.plot(x_ext4[6:23,0], x_ext4[6:23,1],'r', label='Swing Phase')
    plt.title(' Rotational Velocity of Shank vs. Time (Full Swing Phase)')
    plt.xlabel('Time (s)')
    plt.ylabel('Rotational Velocity (rad/s)')
    plt.legend()
    plt.show()
    plt.plot(x_ext3[:, 0], x_ext3[:, 1])
    plt.title('Orientation of Shank vs. Time')
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
