import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from data import Data
import math


def plot_data():
    A = np.asarray(Data.angle_data)
    dt = np.diff(A[:,0])
    dt_m = np.ma.array(dt, mask=(dt==0))
    dy = np.diff(A[:,1])
    dy_m = np.ma.array(dy, mask=(dy==0))
    dydt = dy_m/dt_m
    y = scipy.signal.savgol_filter(dydt, 25, 3)

    d2t = np.diff(dt)
    # print(d2t)
    d2t_m = np.ma.array(d2t, mask=(d2t==0))
    d2y = np.diff(dy)
    # print(d2y)
    d2y_m = np.ma.array(d2y, mask=(d2y==0))
    d2yd2t = d2y_m/d2t_m
    y2 = scipy.signal.savgol_filter(d2yd2t, 51, 3)
    plt.plot(A[12:len(d2t),0], y2[12:])
    plt.show()
    plt.plot(A[12:len(dt),0], y[12:])
    plt.plot(A[12:,0], A[12:,1])
    plt.show()

    radius = 0.4
    


def main():
    plot_data()


if __name__ == '__main__':
    main()
