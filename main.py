import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from data import Data


def plot_data():
    A = np.asarray(Data.angle_data)
    dt = np.diff(A[:,0])
    dt_m = np.ma.array(dt, mask=(dt==0))
    dt_clean = []
    dy = np.diff(A[:,1])
    dy_m = np.ma.array(dy, mask=(dy==0))
    dydt = dy_m/dt_m
    y = scipy.signal.savgol_filter(dydt, 25, 3)
    plt.plot(A[:len(dt),0], y)
    plt.plot(A[:,0], A[:,1])
    plt.show()


def main():
    plot_data()


if __name__ == '__main__':
    main()
