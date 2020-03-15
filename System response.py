# Citation 550 - Linear simulation


from math import *
import numpy as np
import control.matlab as ml
import control
import matplotlib.pyplot as plt
from flight_conditions import *


def plot_params(t, ys) -> None:
    plt.figure(4)
    plt.subplot(221)
    plt.plot(t, ys[0, :])
    plt.title('u rate')
    plt.subplot(222)
    plt.plot(t, ys[1, :])
    plt.title('Alpha rate')
    plt.subplot(223)
    plt.plot(t, ys[2, :])
    plt.title('theta rate')
    plt.subplot(224)
    plt.plot(t, ys[3, :])
    plt.title('q rate')
    plt.tight_layout()
    plt.show()
    return


def plot_response(sys, T=None, U=0.0, impulse=True, step=True) -> None:
    if impulse:
        t, yout = control.impulse_response(sys, T=T, X0=U)
        plot_params(t, yout)
    if step:
        t, yout = control.step_response(sys, T=T, X0=U)
        plot_params(t, yout)
    return


def sym_flight(ac: FlightParams):
    # Construct the symmetric flight state-space system for an aircraft
    C1s = np.matrix([[-2*ac.muc*ac.c/ac.V0, 0 , 0 , 0],
                     [0 , (ac.CZadot -2*ac.muc)*ac.c/ac.V0, 0 , 0],
                     [0 , 0 , -ac.c/ac.V0 , 0],
                     [0, ac.Cmadot , 0 , -2*ac.muc*ac.KY2*ac.c/ac.V0]])

    C2s = np.matrix([[ac.CXu, ac.CXa , ac.CZ0 , ac.CXq],
                     [ac.CZu , ac.CZa, -ac.CX0 , ac.CZq + 2*ac.muc],
                     [0 , 0 , 0 , 1],
                     [ac.Cmu, ac.Cma , 0 , ac.Cmq]])

    C3s = np.matrix([[ac.CXde],
                     [ac.CZde],
                     [0],
                     [ac.Cmde]])

    C1sinv = np.linalg.inv(C1s)

    As = -C1sinv*C2s
    Bs = -C1sinv*C3s
    Cs = np.identity(4)
    Ds =np.zeros([4,1])

    return ml.ss(As , Bs , Cs , Ds)


def asym_flight(ac: FlightParams):
    # Construct the asymmetric flight state-space system for an aircraft
    C1a =np.matrix([[(ac.CYbdot-2*ac.mub), 0 , 0 , 0],
                   [0 , -1/2*ac.b/ac.V0, 0 , 0],
                   [0 , 0 , -4*ac.mub*ac.KX2*ac.b/ac.V0 , 4*ac.mub*ac.KXZ*ac.b/ac.V0],
                   [ac.Cnbdot*ac.b/ac.V0, 0 , 4*ac.mub*ac.KXZ*ac.b/ac.V0 , -4*ac.mub*ac.KX2*ac.b/ac.V0]])

    C2a =np.matrix([[ac.CYb, ac.CL , ac.CYp , ac.CYr-4*ac.mub],
                   [0 , 0, 1 , 0],
                   [ac.Clb , 0 , ac.Clp , ac.Clr],
                   [ac.Cnb, 0 , ac.Cnp , ac.Cnr]])

    C3a = np.matrix([[ac.CYda, ac.CYdr],
                    [0 , 0],
                    [ac.Clda, ac.Cldr],
                    [ac.Cnda, ac.Cndr]])

    C1ainv = np.linalg.inv(C1a)

    Aa = -C1ainv*C2a
    Ba = -C1ainv*C3a
    Ca = np.identity(4)
    Da =np.zeros([4 , 2])

    return ml.ss(Aa , Ba , Ca , Da)


if __name__ == "__main__":
    # Define flying aircraft with default parameters
    aircraft = FlightParams(m=1000)

    # For symmetric flight
    syss = sym_flight(aircraft)
    plot_response(syss, step=True, impulse=True)

    # For asymmetric flight
    sysa = asym_flight(aircraft)
    plot_response(sysa, step=True, impulse=True)

    t = np.arange(0, 10.01, 0.01)
    U1 = np.array([1, 1, 1, 1])
    U0 = np.array(np.zeros((len(t) - len(U1))))
    U = np.concatenate([U1, U0])





