import math
import numpy as np
import csv
import matplotlib.pyplot as plt
import scipy.signal
from data import Data
from ankle_model import FootDropAnkleModel


def write_data(file_name, x_data, y_data):
    with open(file_name, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for i in range(len(x_data)):
            row = [x_data[i],y_data[i]]
            writer.writerow(row)
    # f.close()


def run_simulation():
    """
    Runs simulation for 0.35 s (swing phase). Plots solution (state vector).
    """
    foot_drop = FootDropAnkleModel()
    foot_drop.normalize_times()
    activation_values = [0, 0.25, 0.5, 0.75, 0.95] 
    initial_orientations = [-15, -10, -5]
    duration = 2
    color = ['r', 'b', 'g', 'm', 'c']
    for i in range(len(initial_orientations)):
        soln = foot_drop.simulate(duration, orientation=initial_orientations[i])
        ######## orientation of foot ########
        # data to CSV
        filename_orientation = f"orientation_{initial_orientations[i]}.csv"
        write_data(filename_orientation, (soln.t-1)*0.35, soln.y[1])
        # plot
        plt.plot((soln.t-1)*0.35, soln.y[1], 'b', label="Initial Orientation = %s" % initial_orientations[i])
        plt.title('Orientation of Foot vs. Time')
        plt.ylabel('Orientation (deg)')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show()
        ######## foot velocity rotational ########
        filename_velocity = f"velocity_{initial_orientations[i]}.csv"
        write_data(filename_velocity, (soln.t-1)*0.35, soln.y[2])
        plt.plot((soln.t-1)*0.35, soln.y[2], 'g', label="Initial Orientation = %s" % initial_orientations[i])
        plt.title('Rotational Velocity of Foot vs. Time')
        plt.ylabel('Velocity (deg/s)')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.show()
        print(soln.message)
        print(soln.status)
    for j in range(len(activation_values)):
        ######## activation ########
        soln = foot_drop.simulate(duration, activation=activation_values[j])
        # data to CSV
        filename_activation = f"activation_{activation_values[j]}.csv"
        write_data(filename_activation, (soln.t-1)*0.35, soln.y[0]/2)
        # plot
        plt.subplot(3,2,j+1)
        plt.plot((soln.t-1)*0.35, soln.y[0]/2, color[j], label="Initial Activation Level = % s" % activation_values[j])
        plt.title('Muscle Activation vs. Time')
        plt.legend()
        plt.xlabel('Time (s)')
        plt.ylabel('Activation Level')
    plt.tight_layout()
    plt.show()

    plt.plot((soln.t-1)*0.35, soln.y[0]/2, color[4], label="Initial Activation Level = % s" % activation_values[4])
    plt.title('Muscle Activation vs. Time')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Activation Level')
    plt.show()


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
    # graph_data()


if __name__ == '__main__':
    main()
