# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:38:58 2020

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

#data = loadmat('matlab.mat')
data = loadmat('FTISxprt-20200311_flight3')

FlightParams

c = 2.0569
S = 30
Cm_T_c = -0.0064
d = 27*0.0254
W_s = 60500
rho0 = 1.225

#second stationary measurement series
t = np.array([39*60+16, 40*60+33, 41*60+31, 43*60+15, 44*60+55, 46*60+14, 47*60+8])
hp = np.array([8177, 8432, 8540, 8672, 7835, 7475, 7485])*0.3048
d_e = np.array([0, -0.5, -0.9, -1.45, 0.25, 0.6, 0.85])*pi/180
F_e = np.array([0, -17, -22.5, -41.5, 31.5, 54, 88])
f_used = np.array([708, 728, 742, 772, 792, 810, 826])
T = np.array([2829.76, 2873.19, 2926.55, 2984.99, 2857.65, 2887.34, 2876.39])

#recorded data
t_rec = data['flightdata'][0][0][48][0][0][0].transpose()
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
delta_cg = mcg(867, 2, 2)[2] - mcg(845, 1, 2)[2]
t = 49*60+48
ac.hp = 7408*0.3048
delta_d_e = -0.6*pi/180
f_used = 867

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
p = np.poly1d(np.polyfit(V_e_tilde, d_e_star, 2)) #1 is linear, 2 is parabolic
plt.plot(V_e_tilde,p(V_e_tilde),'r--')
plt.title('Elevator trim curve')
plt.xlabel(r'$\~V_e (m/s)$')
plt.ylabel(r'$\delta^*_e (rad)$')
plt.figure()
plt.plot(V_e_tilde,F_e_star)
q = np.poly1d(np.polyfit(V_e_tilde, F_e_star, 2)) #1 is linear, 2 is parabolic
plt.plot(V_e_tilde,q(V_e_tilde),'r--')
plt.title('Elevator control force curve')
plt.xlabel(r'$\~V_e (m/s)$')
plt.ylabel(r'$F^*_e (N)$')
plt.show()

#calculate Cma
#wat ik oorspronkelijk dacht hoe het moest
#delta_a = data['flightdata'][0][0][16][0][0][0][np.where(t_rec[:,0]==t)]*pi/180
#Cm_alpha = -delta_d_e/delta_a*Cm_d

#formule uit flight dynamics reader
### d_d_e/d_V = 4W/rhoV^3S*1/Cm_d*Cm_alpha/CL_alpha

#met de afgeleide van de trim curve (zit zeker een factor te weinig in)
#deda = np.polyfit(V_e_tilde, d_e_star, 2)[0]
#Cm_alpha = -deda*Cm_d

#d_e vs alpha plotten en daar slope van pakken
alpha = np.array([4, 4, 5, 5, 6, 7, 9])*pi/180
plt.figure()
plt.plot(alpha, d_e_star)
q = np.poly1d(np.polyfit(alpha, d_e_star, 1)) #1 is linear, 2 is parabolic
plt.plot(alpha,q(alpha),'r--')
plt.show()

deda = np.polyfit(alpha, d_e_star, 1)[0]
Cm_alpha = deda*Cm_d