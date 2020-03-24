import numpy as np
from matplotlib import pyplot as plt
import os
from numerical_tools import *


def create_matlabdat(TAT, h, FFL, FFR):
    main_path = os.path.dirname(os.path.realpath(__file__))  # Script folder
    pass


def correct_speed(velocity):
    # Converts from indicated airspeed to calibrated airspeed. Input and output are consistent in units (if input is
    # in knots, output is in knots).
    velocity = velocity * 1.94384449
    tableKIAS = np.append(np.linspace(80,275,40),277)
    tableKCAS = np.array([78, 83, 88, 93, 98, 103, 108, 113, 118, 123, 128, 133, 138, 143, 148, 153, 158, 163, 168, \
                          173, 178, 183, 188, 193, 198, 203, 208, 213, 218, 223, 228, 233, 238, 243, 248, 253, 258, \
                          263, 268, 273, 275])
    multipliers = tableKCAS/tableKIAS
    multipl_func = cont_spline(tableKIAS, multipliers)
    return velocity * multipl_func(velocity).item() * 0.514444444


def correct_mach(mach):
    # Converts from indicated mach number to calibrated mach number.
    table_ind_M = np.append(np.linspace(0.400, 0.700, 31), 0.705)
    table_cal_M = np.append(np.linspace(0.393, 0.693, 31), 0.698)
    multipliers = table_cal_M/table_ind_M
    multipl_func = cont_spline(table_ind_M, multipliers)
    return mach * multipl_func(mach).item()


def ISA_calculator(h, R=287.058, g=9.80665):
    # Takes in an altitude in meters, outputs the temperature, pressure and density in SI units at that altitude.
    h0 = np.array([0, 11000, 20000, 32000, 47000, 51000, 71000, 86000])
    T = np.array([288.15, 216.64999999999998, 216.64999999999998, 228.64999999999998, 270.65, 270.65,
                  214.64999999999998, 184.64999999999998])
    a = np.array([-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002])
    p = np.array([101325, 22625.7915, 5471.93507, 867.254994, 110.766577, 66.8482964, 3.94900089, 0.301619121])
    if h < 0:
        raise ValueError(f"This height is invalid: h = {h} < 0")
    elif h <= h0[1]:
        T1 = T[0] + a[0] * (h - h0[0])
        p1 = p[0] * ((T1 / T[0]) ** (-g / (a[0] * R)))
    elif h <= h0[2]:
        T1 = T[1] + a[1] * (h - h0[1])
        p1 = p[1] * np.exp(-g * (h - h0[1]) / (T1 * R))
    elif h <= h0[3]:
        T1 = T[2] + a[2] * (h - h0[2])
        p1 = p[2] * ((T1 / T[2]) ** (-g / (a[2] * R)))
    elif h <= h0[4]:
        T1 = T[3] + a[3] * (h - h0[3])
        p1 = p[3] * ((T1 / T[3]) ** (-g / (a[3] * R)))
    elif h <= h0[5]:
        T1 = T[4] + a[4] * (h - h0[4])
        p1 = p[4] * np.exp(-g * (h - h0[4]) / (T1 * R))
    elif h <= h0[6]:
        T1 = T[5] + a[5] * (h - h0[5])
        p1 = p[5] * ((T1 / T[5]) ** (-g / (a[5] * R)))
    elif h <= h0[7]:
        T1 = T[6] + a[6] * (h - h0[6])
        p1 = p[6] * ((T1 / T[6]) ** (-g / (a[6] * R)))
    else:
        raise ValueError(f"This height is invalid: {h0[7]} < h = {h}")
    rho = p1 / (R * T1)
    return p1, rho, T1


class AtmosConditions():
    def __init__(self, m = 4000, Ws=60500, hp=1000, V0=200, p0=101325, rho0=1.225, T0=288.15, g=9.80665, R=287.058):
        """Class containing all atmospheric properties of the flight conditions. Changeable parameters are parameters
        ending with a 0. The class automatically updates all other flight parameters. """
        self.gamma = 1.4
        self.g = g
        self.R = R
        self._observers = []
        self._m = m
        self._W = m * g
        self._Ws = Ws
        self._V0 = V0
        self._Vtas = correct_speed(V0)
        self._p0 = p0
        self._rho0 = rho0
        self._T0 = T0
        self._hp = hp
        self._renew_parameters()

    def bind_observer(self, cb):
        self._observers.append(cb)

    def _renew_parameters(self):
        # Updates all relevant flight parameters based on the ISA atmosphere and mach and weight corrections.
        self._p, self._rho, self._T = ISA_calculator(self._hp)
        self._Vtas = correct_speed(self._V0)
        self._M = self._getMach()
        self._T = self._correctT()
        self.a = np.sqrt(self.gamma * self.R * self._T)
        self._V = self._correctV(self._M*self.a)
        for obs in self._observers:
            obs()
        return

    def _getMach(self):
        # Returns the mach number for the given flight conditions.
        return np.sqrt((2 / (self.gamma - 1)) * ((1 + (self._p0 / self._p) * ((1 + ((self.gamma - 1) * self._rho0 * self._Vtas * self._Vtas) /
                                                                               (2 * self.gamma * self._p0)) ** (self.gamma / (self.gamma - 1)) - 1)) ** ((self.gamma - 1) / self.gamma) - 1))

    def _correctT(self):
        # Returns corrected temperature for ram rise.
        return self._T0/(1 + 0.5 * self._M * self._M * (self.gamma - 1))

    def _correctV(self, speed):
        # Returns equivalent airspeed, corrected for atmospheric differences and weight differences.
        return speed*np.sqrt(self._rho/self._rho0)*np.sqrt(self._Ws/self._W)

    @property
    def m(self):
        return self._m

    @property
    def W(self):
        return self._W

    @property
    def Ws(self):
        return self._Ws

    @property
    def V0(self):
        return self._V0

    @property
    def Vtas(self):
        return self._Vtas

    @property
    def p0(self):
        return self._p0

    @property
    def rho0(self):
        return self._rho0

    @property
    def T0(self):
        return self._T0

    @property
    def hp(self):
        return self._hp

    @property
    def V(self):
        return self._V

    @property
    def p(self):
        return self._p

    @property
    def rho(self):
        return self._rho

    @property
    def T(self):
        return self._T

    @property
    def M(self):
        return self._M

    @m.setter
    def m(self, value):
        self._m = value
        self._W = value * self.g
        self._renew_parameters()

    @W.setter
    def W(self, value):
        self._W = value
        self._m = value / self.g
        self._V = self._correctV(self._M*self.a)
        for obs in self._observers:
            obs()

    @Ws.setter
    def Ws(self, value):
        self._Ws=value
        self._renew_parameters()

    @V0.setter
    def V0(self, value):
        self._V0=value
        self._renew_parameters()

    @Vtas.setter
    def Vtas(self, value):
        self._Vtas=value
        self._M = self._getMach()
        self._T = self._correctT()
        self.a = np.sqrt(self.gamma * self.R * self._T)
        self._V = self._correctV(self._M*self.a)
        for obs in self._observers:
            obs()
        return

    @p0.setter
    def p0(self, value):
        self._p0=value
        self._renew_parameters()

    @rho0.setter
    def rho0(self, value):
        self._rho0=value
        self._renew_parameters()

    @T0.setter
    def T0(self, value):
        self._T0=value
        self._renew_parameters()

    @hp.setter
    def hp(self, value):
        self._hp=value
        self._renew_parameters()

    @M.setter
    def M(self, value):
        self._M=correct_mach(value)
        self._T = self._correctT()
        self.a = np.sqrt(self.gamma * self.R * self._T)
        self._V = self._correctV(self._M*self.a)
        for obs in self._observers:
            obs()

    @V.setter
    def V(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was V")

    @p.setter
    def p(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was p")

    @rho.setter
    def rho(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was rho")

    @T.setter
    def T(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was T")


class FlightParams(AtmosConditions):
    def __init__(self, hp=1000, V0=200, alpha0=0, th0=0, m=4000, Ws=60500, e=0.8, CD0=0.04, CLa=5.084, Cma=-0.5626,
                 Cmde=-1.1642, S=30., Sh_multiply=0.2, lh=0.71 * 5.968, c=2.0569, b=15.911, bh=5.791, Vh_V=1,
                 ih=-2 * np.pi / 180, rho0=1.2250, T0=288.15, R=287.058, g=9.80665, KX2=0.019, KZ2=0.042, KXZ=0.002,
                 KY2=1.3925, Cm0=0.0297, CXu=-0.095, CXa=-0.4797, CXadot=0.0833, CXq=-0.2817, CXde=-0.0373,
                 CZu=-0.3762, CZa=-5.7434, CZadot=-0.0035, CZq=-5.6629, CZde=-0.6961, Cmu=0.0699, Cmadot=0.178,
                 Cmq=-8.7941, CmTc=-0.0064, CYb=-0.75, CYbdot=0, CYp=-0.0304, CYr=0.8495, CYda=-0.04, CYdr=0.23,
                 Clb=-0.1026, Clp=-0.7108, Clr=0.2376, Clda=-0.2309, Cldr=0.0344, Cnb=0.1348, Cnbdot=0, Cnp=-0.0602,
                 Cnr=-0.2061, Cnda=-0.012, Cndr=-0.0939):

        # Inherit the parameters gamma, g, R, m, ms, W, Ws, V0, p0, rho0, T0, hp, a, V, p, rho and T from the
        # atmospheric flight conditions class. All 'output' parameters (a, V, p, rho and T) are corrected and reduced.
        super().__init__(m=m, Ws=Ws, hp=hp, V0=V0, rho0=rho0, T0=T0, g=g, R=R)

        self.alpha0 =    alpha0             # angle of attack in the stationary flight condition [rad]
        self.th0    =    th0                # pitch angle in the stationary flight condition [rad]

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

        # Constant values concerning aircraft inertia
        self.KX2    = KX2
        self.KZ2    = KZ2
        self.KXZ    = KXZ
        self.KY2    = KY2


        # Aerodynamic constants
        self.Cmac   = Cm0                                       # Moment coefficient about the aerodynamic centre [ ]
        self.CNwa   = CLa                                       # Wing normal force slope [ ]
        self.CNha   = 2 * np.pi * self.Ah / (self.Ah + 2)       # Stabiliser normal force slope [ ]
        self.depsda = 4 / (self.A + 2)                          # Downwash gradient [ ]

        # Stability derivatives
        self.CXu    = CXu
        self.CXa    = -CXa		# Positive! (has been erroneously negative since 1993)
        self.CXadot = CXadot
        self.CXq    = CXq
        self.CXde   = CXde

        self.CZu    = CZu
        self.CZa    = CZa
        self.CZadot = CZadot
        self.CZq    = CZq
        self.CZde   = CZde

        self.Cmu    = Cmu
        self.Cmadot = Cmadot
        self.Cmq    = Cmq
        self.CmTc   = CmTc

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

        # Make sure parameters are auto-updated
        self.bind_observer(self._renew_stability_param)
        self._renew_stability_param()

    @property
    def muc(self):
        return self._muc

    @property
    def mub(self):
        return self._mub

    @property
    def CL(self):
        return self._CL

    @property
    def CD(self):
        return self._CD

    @property
    def CX0(self):
        return self._CX0

    @property
    def CZ0(self):
        return self._CZ0

    def _renew_stability_param(self):
        self._muc    = self.m / (self.rho * self.S * self.c)
        self._mub    = self.m / (self.rho * self.S * self.b)

        # Lift and drag coefficient
        self._CL_eig = 2 * self.W / (self.rho * self.Vtas * self.Vtas * self.S)               # Lift coefficient [ ]
        self._CL = 2 * self.W / (self.rho * self.V * self.V * self.S)
        self._CD = self.CD0+(self.CLa * self.alpha0)**2 / (np.pi * self.A * self.e) # Drag coefficient [ ]

        self._CX0    = self.W * np.sin(self.th0) / (0.5 * self.rho * self.V * self.V * self.S)
        self._CZ0    = -self.W * np.cos(self.th0) / (0.5 * self.rho * self.V * self.V * self.S)

    @muc.setter
    def muc(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was muc")

    @mub.setter
    def mub(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was mub")

    @CL.setter
    def CL(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was CL")

    @CD.setter
    def CD(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was CD")

    @CX0.setter
    def CX0(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was CX0")

    @CZ0.setter
    def CZ0(self, value):
        print("Changing a corrected parameter directly is not allowed. Change the measured parameter instead.")
        print("Parameter you tried to change was CZ0")







if __name__ == "__main__":
    ac = FlightParams()
    print(ac.V)
    print(ac.CL)
    ac.V0 = 168
    print(ac.V)
    print(ac.CL)
