import numpy as np
from matplotlib import pyplot as plt

class FlightParams():
    def __init__(self, hp0=1000, V0=200, alpha0=0, th0=0, m=1000, e=0.8, CD0=0.04, CLa=5.084, Cma=-0.5626, Cmde=-1.1642,\
                 S=30., Sh_multiply=0.2, lh=0.71*5.968, c=2.0569, b=15.911, bh=5.791, Vh_V=1, ih=-2*np.pi/180, \
                 rho0=1.2250, lamba=-0.0065, Temp0=288.15, R=287.05, g=9.80665, KX2=0.019, KZ2=0.042, KXZ=0.002, \
                 KY2=1.25*1.114, Cmac=0, CXu=-0.02792, CXa=0.47966, CXadot=0.0833, CXq=-0.28170, CXde=-0.03728, \
                 CZu=-0.37616, CZa=-5.74340, CZadot=-0.00350, CZq=-5.66290, CZde=-0.69612, Cmu=0.06990, Cmadot=0.178, \
                 Cmq=-8.79415, CYb=-0.75,CYbdot=0, CYp=-0.0304, CYr=0.8495, CYda=-0.04, CYdr=0.23, Clb=-0.1026, \
                 Clp=-0.71085, Clr=0.2376, Clda=-0.23088, Cldr=0.0344, Cnb=0.1348, Cnbdot=0, Cnp=-0.0602, Cnr=-0.2061, \
                 Cnda=-0.012, Cndr=-0.0939):

        self.hp0    =    hp0                # pressure altitude in the stationary flight condition [m]
        self.V0     =    V0                 # true airspeed in the stationary flight condition [m/sec]
        self.alpha0 =    alpha0             # angle of attack in the stationary flight condition [rad]
        self.th0    =    th0                # pitch angle in the stationary flight condition [rad]

        # Aircraft mass
        self.m      =    m                  # mass [kg]

        # aerodynamic properties
        self.e      =    e                  # Oswald factor [ ]
        self.CD0    =    CD0                # Zero lift drag coefficient [ ]
        self.CLa    =    CLa                # Slope of CL-alpha curve [ ]

        # Longitudinal stability
        self.Cma    =    Cma                # longitudinal stabilty [ ]
        self.Cmde   =    Cmde               # elevator effectiveness [ ]

        # Aircraft geometry
        self.S      = S                     # wing area [m^2]
        self.Sh     = Sh_multiply * S       # stabiliser area [m^2]
        self.Sh_S   = self.Sh / S           # [ ]
        self.lh     = lh                    # tail length [m]
        self.c      = c                     # mean aerodynamic cord [m]
        self.lh_c   = lh / c                # [ ]
        self.b      = b                     # wing span [m]
        self.bh     = bh                    # stabilser span [m]
        self.A      = b * b / S             # wing aspect ratio [ ]
        self.Ah     = bh * bh / self.Sh     # stabilser aspect ratio [ ]
        self.Vh_V   = Vh_V                  # [ ]
        self.ih     = ih                    # stabiliser angle of incidence [rad]

        # Constant values concerning atmosphere and gravity
        self.rho0   = rho0                  # air density at sea level [kg/m^3]
        self.lamba  = lamba                 # temperature gradient in ISA [K/m]
        self.Temp0  = Temp0                 # temperature at sea level in ISA [K]
        self.R      = R                     # specific gas constant [m^2/sec^2K]
        self.g      = g                     # [m/sec^2] (gravity constant)

        # air density [kg/m^3]
        self.rho    = rho0 * pow( ((1+(lamba * hp0 / Temp0))), (-((g / (lamba*R)) + 1)))
        self.W      = m * g                 # [N]       (aircraft weight)

        # Constant values concerning aircraft inertia
        self.muc    = m / (self.rho * S * c)
        self.mub    = m / (self.rho * S * b)
        self.KX2    = KX2
        self.KZ2    = KZ2
        self.KXZ    = KXZ
        self.KY2    = KY2

        # Aerodynamic constants
        self.Cmac   = Cmac                                      # Moment coefficient about the aerodynamic centre [ ]
        self.CNwa   = CLa                                       # Wing normal force slope [ ]
        self.CNha   = 2 * np.pi * self.Ah / (self.Ah + 2)       # Stabiliser normal force slope [ ]
        self.depsda = 4 / (self.A + 2)                          # Downwash gradient [ ]

        # Lift and drag coefficient
        self.CL = 2 * self.W / (self.rho * V0 * V0 * S)         # Lift coefficient [ ]
        self.CD = CD0+(CLa * alpha0)**2 / (np.pi * self.A * e)  # Drag coefficient [ ]

        # Stabiblity derivatives
        self.CX0    = self.W * np.sin(th0) / (0.5 * self.rho * V0 * V0 * S)
        self.CXu    = -CXu
        self.CXa    = CXa		# Positive! (has been erroneously negative since 1993)
        self.CXadot = CXadot
        self.CXq    = CXq
        self.CXde   = CXde

        self.CZ0    = -self.W * np.cos(th0) / (0.5 * self.rho * V0 * V0 * S)
        self.CZu    = CZu
        self.CZa    = CZa
        self.CZadot = CZadot
        self.CZq    = CZq
        self.CZde   = CZde

        self.Cmu    = Cmu
        self.Cmadot = Cmadot
        self.Cmq    = Cmq

        self.CYb    = CYb
        self.CYbdot = CYbdot
        self.CYp    = CYp
        self.CYr    = CYr
        self.CYda   = CYda
        self.CYdr   = CYdr

        self.Clb    = Clb
        self.Clp    = Clp
        self.Clr    = Clr
        self.Clda   = Clda
        self.Cldr   = Cldr

        self.Cnb    =  Cnb
        self.Cnbdot =  Cnbdot
        self.Cnp    =  Cnp
        self.Cnr    =  Cnr
        self.Cnda   =  Cnda
        self.Cndr   =  Cndr