# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:52:56 2020

@author: Gytha
"""

#### Elevator trim curve and elevator control force curve ####

import numpy as np
from numpy import *
from math import *
import scipy
from scipy.io import loadmat
from matplotlib import pyplot as plt
from MCGV2 import *   #choose ref data or flight data
from System_response import *
from flight_conditions import *
from numerical_tools import *

data = loadmat('matlab.mat')
#data = loadmat('FTISxprt-20200311_flight3')

FlightParams

c = 2.0569
S = 30
Cm_T_c = -0.0064
d = 27*0.0254
W_s = 605000

#second stationary measurement series
t = np.array([37*60+19, 39*60+11, 41*60+24, 42*60+56, 45*60+41, 47*60+20, 48*60+40])
hp = np.array([6060, 6350, 6550, 6880, 6160, 5810, 5310])*0.3048
d_e = np.array([0, -0.4, -0.9, -1.5, 0.4, 0.6, 1])*pi/180
F_e = np.array([0, -23, -29, -46, 26, 40, 83])
f_used = np.array([664, 694, 730, 755, 798, 825, 846])
T = np.array([3676.79, 3688.19, 3709.39, 3712.83, 3553.91, 3670.39, 3793.06])

#recorded data
t_rec = data['flightdata'][0][0][47][0][0][0].transpose()
V = Dadc1_tas = data['flightdata'][0][0][41][0][0][0]*0.514444444

W = np.zeros(7)
V_t = np.zeros(7)
V_e_tilde = np.zeros(7)
F_e_star = np.zeros(7)
for i in range(7):
    ac.hp = hp[i]
    
    #compute weight
    W[i] = mcg(f_used[i], 1, 2)[0]
    
    #compute reduced airspeed
    V_t[i] = V[np.where(t_rec[:,0]==t[i])]
    V_e_tilde[i] = V_t[i]*np.sqrt(ac.rho/rho0)*np.sqrt(W_s/W[i])
    
    #compute reduced elevator control force
    F_e_star[i] = F_e[i]*W_s/W[i]

#gravity shift -> new parameters
delta_cg = mcg(910, 2, 2)[2] - mcg(881, 0, 2)[2]
t = 52*60+46
ac.hp = 5790*0.3048
delta_d_e = -0.5*pi/180
f_used = 910

W = mcg(f_used, 2, 2)[0]
CL = 2*W/(ac.rho*(V[np.where(t_rec[:,0]==t)])**2*S)

#compute reduced elevator deflection
Cm_d = -1/delta_d_e*CL*delta_cg/c
Tc_s = np.zeros(7)
Tc = np.zeros(7)
d_e_star = np.zeros(7)
for i in range(7):
    ac.hp = hp[i]
    Tc[i] = 2*T[i]/(ac.rho*V_t[i]**2*S)
    Tc_s[i] = 2*T[i]/(ac.rho*V_t[i]**2*d**2)
    d_e_star[i] = d_e[i]-Cm_T_c/Cm_d*(Tc_s[i]-Tc[i])
    
#sort data
V_e_tilde = np.sort(V_e_tilde)
d_e_star = np.sort(d_e_star)
F_e_star = np.sort(F_e_star)

# plot curves
plt.figure()
plt.plot(V_e_tilde,d_e_star)
p = np.poly1d(np.polyfit(V_e_tilde, d_e_star, 2))
plt.plot(V_e_tilde,p(V_e_tilde),'r--')
plt.title('Elevator trim curve')
plt.xlabel(r'$\~V_e$')
plt.ylabel(r'$\delta^*_e$')
plt.figure()
plt.plot(V_e_tilde,F_e_star)
q = np.poly1d(np.polyfit(V_e_tilde, F_e_star, 2))
plt.plot(V_e_tilde,q(V_e_tilde),'r--')
plt.title('Elevator control force curve')
plt.xlabel(r'$\~V_e$')
plt.ylabel(r'$F^*_e$')
plt.show()

#calculate Cma
alpha = np.array([3.4, 4.1, 4.5, 5.3, 6.3, 7.3, 8.5])*pi/180
plt.figure()
plt.plot(alpha, d_e_star)
q = np.poly1d(np.polyfit(alpha, d_e_star, 1)) #1 is linear, 2 is parabolic
plt.plot(alpha,q(alpha),'r--')
plt.show()

deda = np.polyfit(alpha, d_e_star, 1)[0]
Cm_alpha = deda*Cm_d