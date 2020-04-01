import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from data import Data
from ankle_model import FootDropAnkleModel
import math


def plot_data():
    #### NOT BEING USED CURRENTLY - KIRA APR 1
    A = np.asarray(Data.angle_data)
    angle_rad = []
    for val in A[:,1]:
        angle_rad.append(val*(math.pi/180))

    dt = np.diff(A[:,0])
    dt_m = np.ma.array(dt, mask=(dt==0))
    dy = np.diff(angle_rad)
    dy_m = np.ma.array(dy, mask=(dy==0))
    dydt = dy_m/dt_m
    y = scipy.signal.savgol_filter(dydt, 25, 3)

    d2t = np.diff(dt)
    d2t_m = np.ma.array(d2t, mask=(d2t==0))
    d2y = np.diff(dy)
    d2y_m = np.ma.array(d2y, mask=(d2y==0))
    d2yd2t = d2y_m/d2t_m
    y2 = scipy.signal.savgol_filter(d2yd2t, 51, 3)
    div_data = []
    # for value in y2:
    #     div_data.append(value/10000000)
    plt.plot(A[12:len(d2t),0], y2[12:])
    # plt.show()
    plt.plot(A[12:len(dt),0], y[12:])

    # plt.plot(A[12:len(d2t),0], div_data[12:])
    plt.plot(A[12:,0], A[12:,1])
    plt.show()

    radius = 0.4


def run_simulation():
    foot_drop = FootDropAnkleModel()
    # foot_drop.simulate(4.0)
    # foot_drop.simulate_rk4()
    soln = foot_drop.simulate(0.35)
    plt.plot(soln.t, soln.y[0], 'r')
    plt.show()
    # print(foot_drop.muscle_ex_data)
    plt.plot(soln.t, soln.y[1], 'b')
    plt.show()
    plt.plot(soln.t, soln.y[2], 'g')
    plt.show() 
    # print(soln.t)
    print(soln.message)
    print(soln.status)

def graph_data():
    x_ext1 = np.array(Data.linear_acc_shank_x)
    x_ext2 = np.array(Data.linear_acc_shank_z)
    x_ext3 = np.array(Data.abs_orientation_shank) 
    x_ext4 = np.array(Data.abs_velocity_rotation_shank)
    u1 = np.array(Data.muscle_excitation_level_fig3)
    plt.plot(u1[:,0], u1[:,1])
    # plt.plot(x_ext1[:,0], x_ext1[:,1], label='x-direction')
    # plt.plot(x_ext2[:,0], x_ext2[:,1], label='z-direction')
    # plt.title('Linear Acceleration of the Shank vs. Time')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Acceleration (m/s^2)')
    # plt.show()
    # plt.plot(x_ext4[:,0], x_ext4[:,1])
    # plt.title('Absolute Rotational Velocity of Shank vs. Time')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Rotational Velocity (rad/s)')
    # plt.show()
    # plt.plot(x_ext3[:,0], x_ext3[:,1] *(math.pi/180))
    # plt.title('Absolute Orientation of Shank vs. Time')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Orientation (rad)')
    plt.show()

def main():
    run_simulation()
    # graph_data()


if __name__ == '__main__':
    main()
