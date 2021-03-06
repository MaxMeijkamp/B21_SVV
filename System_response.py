# Citation 550 - Linear simulation


from math import *
import numpy as np
import control.matlab as ml
import control
import matplotlib.pyplot as plt
from flight_conditions import *
from Cit_par_ref import *
from Cit_par_B21 import *
from Cit_par_final import *


def plot_params(t, ys, sym=True) -> None:
    if sym:
        labels = ['u', 'alpha', 'theta', 'q']
    else:
        labels = ['Sideslip angle', 'Roll angle', 'Roll rate', 'Yaw rate']
    plt.figure(4)
    plt.subplot(221)
    plt.plot(t, ys[0, :])
    plt.title(labels[0])
    plt.subplot(222)
    plt.plot(t, ys[1, :])
    plt.title(labels[1])
    plt.subplot(223)
    plt.plot(t, ys[2, :])
    plt.title(labels[2])
    plt.subplot(224)
    plt.plot(t, ys[3, :])
    plt.title(labels[3])
    plt.tight_layout()
    plt.show()
    return


def plot_response(sys, T=None, U=0.0, X0=np.zeros(4), impulse=False, step=False, lsim=False, sym=True, input=0) -> None:
    if impulse:
        t, yout = control.impulse_response(sys, T=T, X0=X0, input=0)
        plot_params(t, yout, sym=sym)
    if step:
        t, yout = control.step_response(sys, T=T, X0=X0, input=0)
        plot_params(t, yout, sym=sym)
    if lsim:
        t, yout, xout = control.forced_response(sys, U=U, T=T, X0=X0)
        plot_params(t, yout, sym=sym)
    return


def inspect_sys(sys, eigvals=True, eigvect=True) -> None:
    vals, vects = np.linalg.eig(sys.A)
    if eigvals:
        print("First pair = {}, norm = {}".format(vals[0:2], np.linalg.norm(vals[0])), sep='')
        print("Second pair = {}, norm = {}".format(vals[2:], np.linalg.norm(vals[2])), sep='')
    if eigvect:
        print("First & second vect = {} and {}".format(vects[:,0], vects[:,1]))
        print("Third & fourth vect = {} and {}".format(vects[:,2], vects[:,3]))


def sym_flight(ac: FlightParams):
    # Construct the symmetric flight state-space system for an aircraft
    C1s = np.matrix([[-2*ac.muc*ac.c/ac.Vtas, 0 , 0 , 0],
                     [0 , (ac.CZadot -2*ac.muc)*ac.c/ac.Vtas, 0 , 0],
                     [0 , 0 , -ac.c/ac.Vtas , 0],
                     [0, ac.Cmadot*ac.c/ac.Vtas , 0 , -2*ac.muc*ac.KY2*ac.c/ac.Vtas]])

    C2s = np.matrix([[-ac.CXu, -ac.CXa , -ac.CZ0 , -ac.CXq],
                     [-ac.CZu , -ac.CZa, ac.CX0 , -ac.CZq - 2*ac.muc],
                     [0 , 0 , 0 , -1],
                     [-ac.Cmu, -ac.Cma , 0 , -ac.Cmq]])

    C3s = np.matrix([[-ac.CXde],
                     [-ac.CZde],
                     [0],
                     [-ac.Cmde]])

    C1sinv = np.linalg.inv(C1s)

    As = C1sinv*C2s
    Bs = C1sinv*C3s
    Cs = np.identity(4)
    Ds =np.zeros([4,1])

    return ml.ss(As, Bs, Cs, Ds)


def asym_flight(ac: FlightParams):
    # Construct the asymmetric flight state-space system for an aircraft
    C1a =np.matrix([[(ac.CYbdot-2*ac.mub)*ac.b/ac.V, 0 , 0 , 0],
                   [0 , -0.5*ac.b/ac.V, 0 , 0],
                   [0 , 0 , -4*ac.mub*ac.KX2*ac.b/ac.V , 4*ac.mub*ac.KXZ*ac.b/ac.V],
                   [ac.Cnbdot*ac.b/ac.V, 0 , 4*ac.mub*ac.KXZ*ac.b/ac.V , -4*ac.mub*ac.KZ2*ac.b/ac.V]])

    C2a =np.matrix([[-ac.CYb, -ac.CL , -ac.CYp , -ac.CYr+4*ac.mub],
                   [0 , 0, -1 , 0],
                   [-ac.Clb , 0 , -ac.Clp , -ac.Clr],
                   [-ac.Cnb, 0 , -ac.Cnp , -ac.Cnr]])

    C3a = np.matrix([[-ac.CYda, -ac.CYdr],
                    [0 , 0],
                    [-ac.Clda, -ac.Cldr],
                    [-ac.Cnda, -ac.Cndr]])

    C1ainv = np.linalg.inv(C1a)

    Aa = C1ainv*C2a
    Ba = C1ainv*C3a
    Ca = np.identity(4)
    Ca[1,:] = Aa[1,:]
    Da = np.zeros([4 , 2])
    Da[1,:] = Ba[1,:]
    return ml.ss(Aa , Ba , Ca , Da)


if __name__ == "__main__":
    # Define flying aircraft with default parameters
    from MCG import *
    W, _, _ = mcg(946.8, 0, 1)
    aircraft = ac_ref(m=W/9.80665)
    aircraft2 = ac_B21(m=W/9.80665)

    # Define input and time-line for reponse
    t = np.arange(0, 100.01, 0.01)
    U = np.zeros((2,t.size))
    U[0,np.where(t<5.)] = np.zeros_like(np.where(t<5.))
    # U1[:,np.where(t1<0.5)] = 1.0

    # For symmetric flight
    syssref = sym_flight(aircraft)
    syssb21 = sym_flight(aircraft2)
    # t, yout = control.impulse_response(syss)
    # plt.plot(t, yout[0, :])
    # plt.show()
    # inspect_sys(syss)
    # plot_response(syss, impulse=True)
    # plot_response(syss, step=True, T=t, U=0.0)
    # plot_response(syssref, lsim=True, T=t, X0=np.zeros(4), U=U[1,:])
    plot_response(syssref, lsim=True, T=t, X0=np.array([0.1, 0, 0, 0]))
    # plot_response(syssb21, lsim=True, T=t, X0=np.zeros(4), U=U[1, :])
    # plot_response(syssb21, step=True, impulse=True, T=t, X0=np.zeros(4))

    # For asymmetric flight
    # sysa = asym_flight(aircraft)
    # plot_response(sysa, impulse=True, T=t, U=U, sym=False)
    # plot_response(sysa, lsim=True, T=t, X0=np.zeros(4), U=U, sym=False)





