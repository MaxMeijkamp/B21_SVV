# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 09:46:36 2020

@author: Gytha
"""

from numpy import *
from flight_conditions import *

# Citation 550 - Linear simulation

# xcg = 0.25 * c

# Stationary flight condition
# taken right before phugoid
def ac_B21(m=6000) -> FlightParams:
    hp0    =  7567*0.3048     	      # pressure altitude in the stationary flight condition [m]
    V0     = 151*0.51444444444            # true airspeed in the stationary flight condition [m/sec]
    alpha0 =  6*pi/180           # angle of attack in the stationary flight condition [rad]
    th0    = 10*pi/180            # pitch angle in the stationary flight condition [rad]
    gamma0 = th0 - alpha0

    # Aircraft mass

    m_cur  =   m          # mass [kg]

    # aerodynamic properties #NOG FOUT
    e      = 0.8          # Oswald factor [ ]
    CD0    = 0.0182898           # Zero lift drag coefficient [ ]
    CLa    =  0.12523*180/np.pi        # Slope of CL-alpha curve [ ]

    # Longitudinal stability
    Cma    =  -0.573          # longitudinal stabilty [ ]
    Cmde   =   -1.138          # elevator effectiveness [ ]

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
    labda = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81            # [m/sec^2] (gravity constant)

    # air density [kg/m^3]
    rho    = rho0 * power( ((1+(labda * hp0 / Temp0))), (-((g / (labda*R)) + 1)))
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
    CXu    = 2 * CL * tan(gamma0)
    CXa    = CL * (1 - 2 * CLa / (pi * A * e))		# Positive! (has been erroneously negative since 1993)
    CXadot = 0
    CXq    = 0
    CXde   = 0

    CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = - 2 * CL
    CZa    = - CLa
    CZadot = - CNha * Vh_V ** 2 * Sh_S * lh_c * depsda
    CZq    = - 2 * CNha * Vh_V ** 2 * Sh_S * lh_c
    CZde   = -0.69612

    Cmu    = 0
    Cmadot = - CNha * Vh_V ** 2 * Sh_S * lh_c ** 2 * depsda
    Cmq    = - 1.1 * CNha * Vh_V ** 2 * Sh_S * lh_c

    CYb    = - 0.46
    CYbdot =  0
    CYp    = - 0.91
    CYr    = - 0.15
    CYda   = 0
    CYdr   = 0.2300

    Clb    = - 0.11
    Clp    = - 0.55
    Clr    = 0.088
    Clda   = -0.23088
    Cldr   = 0.03440

    Cnb    =  0.086
    Cnbdot =   0
    Cnp    = - 0.0097
    Cnr    = - 0.15
    Cnda   = -0.0120
    Cndr   =  -0.0939
    return FlightParams(hp=hp0, V0=V0, alpha0=alpha0, th0=th0, m=m_cur, CD0=CD0, CLa=CLa, Cm0=Cmac, CXu=CXu, CXa=CXa,
                        CXadot=CXadot, CXq=CXq, CXde=CXde, CZu=CZu, CZa=CZa, CZadot=CZadot, CZq=CZq, CZde=CZde, Cmu=Cmu,
                        Cmadot=Cmadot, Cmq=Cmq, CYb=CYb, CYbdot=CYbdot, CYp=CYp, CYr=CYr, CYda=CYda, CYdr=CYdr, Clb=Clb,
                        Clp=Clp, Clr=Clr, Clda=Clda, Cldr=Cldr, Cnb=Cnb, Cnbdot=Cnbdot, Cnp=Cnp, Cnr=Cnr, Cnda=Cnda,
                        Cndr=Cndr, Cma=Cma, Cmde=Cmde)

if __name__ == "__main__":
    from MCG import *
    W, _, _ = mcg(946.8, 0, 1)
    ac = ac_B21(m=W/9.80665)
    ac1 = FlightParams(m=ac.m, hp=ac.hp, V0=ac.V0, alpha0=ac.alpha0, th0=ac.th0)
    print_acdata(ac)
