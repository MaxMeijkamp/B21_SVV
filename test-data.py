# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:14:11 2020

Author: Sander Orbons
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt

dat = loadmat('FTISxprt-20200309_flight4.mat')

# Check Data
#for i in range(49):
#    print(i)
#    print(dat['flightdata'][0][0][i][0])
#   print('_________________')

# Asymmetric Responses
    # for Dutch Roll, Aperiodic Roll and Spiral



def Asym(t_start,t_end):
    
    # Making empty lists
    t_list = []
    pitchr_list = []
    rollr_list = []
    yawr_list = []
    ur_list = []
    qr_list = []
    rr_list = []
    # Appending data for each datapoint
    for i in range(10*t_start,10*t_end):
        # Reading data from file
        t = dat['flightdata'][0][0][48][0][0][0][0][i]
        pitchr = dat['flightdata'][0][0][28][0][0][0][i][0]
        rollr = dat['flightdata'][0][0][27][0][0][0][i][0]
        yawr = dat['flightdata'][0][0][29][0][0][0][i][0]
        ur = dat['flightdata'][0][0][30][0][0][0][i][0]
        qr = dat['flightdata'][0][0][31][0][0][0][i][0]
        rr = dat['flightdata'][0][0][32][0][0][0][i][0]
        # Appending data
        t_list.append(t)
        pitchr_list.append(pitchr)
        rollr_list.append(rollr)
        yawr_list.append(yawr)
        ur_list.append(ur)
        qr_list.append(qr)
        rr_list.append(rr)
    
    # Plotting figures for angle response
    f1 = plt.subplot(3, 1, 1)
    plt.plot(t_list, pitchr_list)
    plt.setp(f1.get_xticklabels(), visible=False)
    plt.title('Pitch, Roll and Yaw Response')
    plt.ylabel('Pitch Rate')
    
    f2 = plt.subplot(3, 1, 2)
    plt.plot(t_list, rollr_list)
    plt.setp(f2.get_xticklabels(), visible=False)
    plt.ylabel('Roll Rate')
    
    f3 = plt.subplot(3, 1, 3)
    plt.plot(t_list, yawr_list)
    plt.setp(f3.get_xticklabels(), fontsize=6)
    plt.ylabel('Yaw Rate')
    plt.xlabel('Time [s]')
    plt.show()
    
    # Plotting figures for  response
    f1 = plt.subplot(3, 1, 1)
    plt.plot(t_list, ur_list)
    plt.setp(f1.get_xticklabels(), visible=False)
    plt.title('u,q and r Response')
    plt.ylabel('u Rate')
    
    f2 = plt.subplot(3, 1, 2)
    plt.plot(t_list, qr_list)
    plt.setp(f2.get_xticklabels(), visible=False)
    plt.ylabel('q Rate')
    
    f3 = plt.subplot(3, 1, 3)
    plt.plot(t_list, rr_list)
    plt.setp(f3.get_xticklabels(), fontsize=6)
    plt.xlabel('Time [s]')
    plt.ylabel('r Rate')
    plt.show()
    return 

#____________________________________________________________________________#
    
# Symmetric response
def Sym(t_start,t_end):
    
    # Making empty lists
    t_list = []
    pitchr_list = []
    TAS_list = []
    ur_list = []
    qr_list = []
    elevator_list = []
    # Appending data for each datapoint
    for i in range(10*t_start,10*t_end):
        # Reading data from file
        t = dat['flightdata'][0][0][48][0][0][0][0][i]
        pitchr = dat['flightdata'][0][0][28][0][0][0][i][0]
        TAS = dat['flightdata'][0][0][42][0][0][0][i][0]
        ur = dat['flightdata'][0][0][30][0][0][0][i][0]
        qr = dat['flightdata'][0][0][31][0][0][0][i][0]
        elevator = dat['flightdata'][0][0][18][0][0][0][i][0]
        # Appending data
        t_list.append(t)
        pitchr_list.append(pitchr)
        TAS_list.append(TAS)
        ur_list.append(ur)
        qr_list.append(qr)
        elevator_list.append(elevator)
        
        
    # Plotting figures for angle response
    f0 = plt.subplot(5, 1, 1)
    plt.plot(t_list, elevator_list)
    plt.setp(f0.get_xticklabels(), visible=False)
    plt.title('Elevator Input')
    plt.ylabel('Deflection Elevator')
    plt.show()
    
    f1 = plt.subplot(5, 1, 2)
    plt.plot(t_list, pitchr_list)
    plt.setp(f1.get_xticklabels(), visible=False)
    plt.title('Response')
    plt.ylabel('Pitch Rate')
    
    f2 = plt.subplot(5, 1, 3)
    plt.plot(t_list, TAS_list)
    plt.setp(f2.get_xticklabels(), visible=False)
    plt.ylabel('TAS')
    plt.show()
    
    f3 = plt.subplot(5, 1, 4)
    plt.plot(t_list, ur_list)
    plt.setp(f3.get_xticklabels(), visible=False)
    plt.ylabel('u Rate')
    
    f4 = plt.subplot(5, 1, 5)
    plt.plot(t_list, qr_list)
    plt.setp(f4.get_xticklabels(), fontsize=6)
    plt.xlabel('Time [s]')
    plt.ylabel('q Rate')
    return
# Phugoid
#print('PHUGOID')
#Sym(55*60,57*60)
## Short Period
#print('SHORT PERIOD')
#Sym(52*60,53*60)
## Dutch Roll
#Asym(50*60 +40, 51*60)
## Dutch Roll YD
#Asym(51*60 +13, 52*60)
## Aperiodic Roll
#Asym(53*60 +30,54*60)
# Spiral
#Asym(60*60,65*60)

