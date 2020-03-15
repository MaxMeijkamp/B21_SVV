# Citation 550 - Linear simulation


from math import *
import numpy as np
import control.matlab as ml
import control
import matplotlib.pyplot as plt
 #xcg = 0.25 * c

# Stationary flight condition

hp0    =    1000   	      # pressure altitude in the stationary flight condition [m]
V0     =   200          # true airspeed in the stationary flight condition [m/sec]
alpha0 =    0         # angle of attack in the stationary flight condition [rad]
th0    =    0         # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =   1000          # mass [kg]

# aerodynamic properties
e      =     0.8        # Oswald factor [ ]
CD0    =   0.04          # Zero lift drag coefficient [ ]
CLa    =      5.084       # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    =  -0.5626          # longitudinal stabilty [ ]
Cmde   =   -1.1642          # elevator effectiveness [ ]

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lamba = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
rho    = rho0 * pow( ((1+(lamba * hp0 / Temp0))), (-((g / (lamba*R)) + 1)))   
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993) 
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939





#For symmetric flight

#t = np.arange(0, 10.01, 0.01)
U1 =np.array([1,1,1,1])
U0 = np.array(np.zeros((len(t)-len(U1))))
U = np.concatenate([U1,U0])
def sym_flight(t,U):
    
    #with t timestep and U input array
    V=V0
    C1s = np.matrix([[-2*muc*c/V, 0 , 0 , 0],
                     [0 , (CZadot -2*muc)*c/V, 0 , 0],
                     [0 , 0 , -c/V , 0],
                     [0, Cmadot , 0 , -2*muc*KY2*c/V]])

    C2s = np.matrix([[CXu, CXa , CZ0 , CXq],
                     [CZu , CZa, -CX0 , CZq + 2*muc],
                     [0 , 0 , 0 , 1],
                     [Cmu, Cma , 0 , Cmq]])

    C3s = np.matrix([[CXde],
                     [CZde],
                     [0],
                     [Cmde]])


    C1sinv = np.linalg.inv(C1s)

    As = -C1sinv*C2s
    Bs = -C1sinv*C3s
    Cs = np.identity(4)
    Ds =np.zeros([4,1])



    syss = ml.ss(As , Bs , Cs , Ds)



    ys,t = control.impulse_response(syss)

    plt.figure(4)
    plt.subplot(221)
    plt.plot(t,ys[:,0])
    plt.title('u rate')
    plt.subplot(222)
    plt.plot(t,ys[:,1])
    plt.title('Alpha rate')
    plt.subplot(223)
    plt.plot(t,ys[:,2])
    plt.title('theta rate')
    plt.subplot(224)
    plt.plot(t,ys[:,3])
    plt.title('q rate')
    plt.show()
    
    return

#For asymmetric flight

V=V0
C1a =np.matrix([[(CYbdot-2*mub), 0 , 0 , 0],
               [0 , -1/2*b/V, 0 , 0],
               [0 , 0 , -4*mub*KX2*b/V , 4*mub*KXZ*b/V],
               [Cnbdot*b/V, 0 , 4*mub*KXZ*b/V , -4*mub*KX2*b/V]])

C2a =np.matrix([[CYb, CL , CYp , CYr-4*mub],
               [0 , 0, 1 , 0],
               [Clb , 0 , Clp , Clr],
               [Cnb, 0 , Cnp , Cnr]])

C3a = np.matrix([[CYda, CYdr],
                [0 , 0],
                [Clda, Cldr],
                [Cnda, Cndr]])



C1ainv = np.linalg.inv(C1a)

Aa = -C1ainv*C2a
Ba = -C1ainv*C3a
Ca = np.identity(4)
Da =np.zeros([4 , 2])


sysa = ml.ss(Aa , Ba , Ca , Da)







