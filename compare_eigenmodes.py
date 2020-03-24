# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:14:11 2020

Author: Sander Orbons
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np
from System_response import *
from MCG import *
from Cit_par_B21 import *


def get_ftis_data(filename):
    data = loadmat(filename)['flightdata']

    keys = np.array(['vane_AOA', 'elevator_dte', 'column_fe', 'lh_engine_FMF', 'rh_engine_FMF', 'lh_engine_itt',
                     'rh_engine_itt', 'lh_engine_OP', 'rh_engine_OP', 'column_Se', 'lh_engine_fan_N1',
                     'lh_engine_turbine_N2', 'rh_engine_fan_N1', 'rh_engine_turbine_N2', 'lh_engine_FU', 'rh_engine_FU',
                     'delta_a', 'delta_e', 'delta_r', 'Gps_date', 'Gps_utcSec', 'Ahrs1_Roll', 'Ahrs1_Pitch',
                     'Fms1_trueHeading', 'Gps_lat', 'Gps_long', 'Ahrs1_bRollRate', 'Ahrs1_bPitchRate', 'Ahrs1_bYawRate',
                     'Ahrs1_bLongAcc', 'Ahrs1_bLatAcc', 'Ahrs1_bNormAcc', 'Ahrs1_aHdgAcc', 'Ahrs1_xHdgAcc', 'Ahrs1_VertAcc',
                     'Dadc1_sat', 'Dadc1_tat', 'Dadc1_alt', 'Dadc1_bcAlt', 'Dadc1_bcAltMb', 'Dadc1_mach', 'Dadc1_cas',
                     'Dadc1_tas', 'Dadc1_altRate', 'measurement_running', 'measurement_n_rdy', 'display_graph_state',
                     'display_active_screen', 'time'])

    ftis = dict.fromkeys(keys)
    units = dict.fromkeys(keys)

    for i in range(keys.size):
        ftis[keys[i]] = data[0,0][i][0,0][0].flatten()
        pre_unit = data[0,0][i][0,0][1]
        while not isinstance(pre_unit, str):
            if np.array(pre_unit).size == 0:
                pre_unit = "-"
                break
            pre_unit = pre_unit[0]
        units[keys[i]] = pre_unit
    return ftis, units

# Check Data
#for i in range(49):
#    print(i)
#    print(dat['flightdata'][0][0][i][0])
#   print('_________________')

# Asymmetric Responses
    # for Dutch Roll, Aperiodic Roll and Spiral


def plot_response(t, data1, data2, data3, data4, label1, label2, label3, label4):
    f0 = plt.subplot(2, 2, 1)
    plt.plot(t, data1)
    plt.setp(f0.get_xticklabels(), visible=False)
    plt.ylabel(label1)

    f1 = plt.subplot(2, 2, 2)
    plt.plot(t, data2)
    plt.setp(f1.get_xticklabels(), visible=False)
    plt.ylabel(label2)

    f2 = plt.subplot(2, 2, 3)
    plt.plot(t, data3)
    plt.setp(f2.get_xticklabels(), visible=False)
    plt.ylabel(label3)

    f3 = plt.subplot(2, 2, 4)
    plt.plot(t, data4)
    plt.setp(f3.get_xticklabels(), fontsize=6)
    plt.ylabel(label4)
    plt.xlabel('Time [s]')
    return


def Asym_2(t_start, t_end, ftis):

    # Gathering data
    t = ftis['time']
    pitchr_list = ftis['Ahrs1_bPitchRate'][(t>t_start) & (t<t_end)]
    rollr_list = ftis['Ahrs1_bRollRate'][(t>t_start) & (t<t_end)]
    yawr_list = ftis['Ahrs1_bYawRate'][(t>t_start) & (t<t_end)]
    ur_list = ftis['Ahrs1_bLongAcc'][(t>t_start) & (t<t_end)]
    qr_list = ftis['Ahrs1_bLatAcc'][(t>t_start) & (t<t_end)]
    rr_list = ftis['Ahrs1_bNormAcc'][(t>t_start) & (t<t_end)]
    t_list = t[(t>t_start) & (t<t_end)]

    # Plotting figures for angle response
    # plot_response(t_list, pitchr_list, rollr_list, , pitchr_list, "TAS [m/s]", "AOA [deg]", "Pitch angle [deg]", "Pitch rate [deg/s]")


def Asym_1(t_start, t_end, ftis):

    # Gathering data
    t = ftis['time']
    pitchr_list = ftis['Ahrs1_bPitchRate'][(t>t_start) & (t<t_end)]
    rollr_list = ftis['Ahrs1_bRollRate'][(t>t_start) & (t<t_end)]
    yawr_list = ftis['Ahrs1_bYawRate'][(t>t_start) & (t<t_end)]
    ur_list = ftis['Ahrs1_bLongAcc'][(t>t_start) & (t<t_end)]
    qr_list = ftis['Ahrs1_bLatAcc'][(t>t_start) & (t<t_end)]
    rr_list = ftis['Ahrs1_bNormAcc'][(t>t_start) & (t<t_end)]
    t_list = t[(t>t_start) & (t<t_end)]
    
    # Plotting figures for angle response
    fig, ax = plt.subplots()
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

    ax.grid(True, which='both')
    # plt.tight_layout()
    plt.show()
    return 

#____________________________________________________________________________#
    
# Symmetric response


def Sym_1(t_start,t_end, ftis):

    t = ftis['time']
    pitchr_list = ftis['Ahrs1_bPitchRate'][(t > t_start) & (t < t_end)]
    TAS_list = ftis['Dadc1_tas'][(t>t_start) & (t<t_end)]
    ur_list = ftis['Ahrs1_bLongAcc'][(t>t_start) & (t<t_end)]
    qr_list = ftis['Ahrs1_bLatAcc'][(t>t_start) & (t<t_end)]
    elevator_list = ftis['delta_e'][(t>t_start) & (t<t_end)]
    t_list = t[(t>t_start) & (t<t_end)]

    # Plotting figures for angle response
    fig, ax = plt.subplots()
    f0 = plt.subplot(5, 1, 1)
    plt.plot(t_list, elevator_list)
    plt.setp(f0.get_xticklabels(), visible=False)
    plt.title('Elevator Input')
    plt.ylabel('Deflection Elevator')
    
    f1 = plt.subplot(5, 1, 2)
    plt.plot(t_list, pitchr_list)
    plt.setp(f1.get_xticklabels(), visible=False)
    plt.title('Response')
    plt.ylabel('Pitch Rate')
    
    f2 = plt.subplot(5, 1, 3)
    plt.plot(t_list, TAS_list)
    plt.setp(f2.get_xticklabels(), visible=False)
    plt.ylabel('TAS')
    
    f3 = plt.subplot(5, 1, 4)
    plt.plot(t_list, ur_list)
    plt.setp(f3.get_xticklabels(), visible=False)
    plt.ylabel('u Rate')
    
    f4 = plt.subplot(5, 1, 5)
    plt.plot(t_list, qr_list)
    plt.setp(f4.get_xticklabels(), fontsize=6)
    plt.xlabel('Time [s]')
    plt.ylabel('q Rate')

    ax.grid(True, which='both')
    # plt.tight_layout()
    plt.show()
    return


def Sym_2(t_start,t_end, ftis, ac: FlightParams):

    t = ftis['time']
    aoa_list = ftis['vane_AOA'][(t>=t_start) & (t<t_end)]
    pitch_list = ftis['Ahrs1_Pitch'][(t>=t_start) & (t<t_end)]
    pitchr_list = ftis['Ahrs1_bPitchRate'][(t >= t_start) & (t < t_end)]
    TAS_list = ftis['Dadc1_tas'][(t>=t_start) & (t<t_end)]
    pitchr_list = np.deg2rad(pitchr_list) * ac.c / TAS_list[0] # Reduce and make dimensionless
    aoa_list = np.deg2rad(aoa_list - aoa_list[0]) # Reduce and change units
    pitch_list = np.deg2rad(pitch_list - pitch_list[0]) # Reduce and change units
    TAS_list = (TAS_list-TAS_list[0])/TAS_list[0] # Make dimensionless and reduce
    t_list = t[(t>=t_start) & (t<t_end)]

    # Plotting figures for angle response

    plot_response(t_list, TAS_list, aoa_list, pitch_list, pitchr_list, "TAS [m/s]", "AOA [deg]", "Pitch angle [deg]", "Pitch rate [deg/s]")
    return


if __name__ == "__main__":
    ftis, units = get_ftis_data('FTISxprt-20200311_flight3.mat')
    # ftis, units = get_ftis_data('matlab.mat')
    tstart = 55*60 # [s]
    tend = 57*60 # [s]
    fuel_used = ftis['lh_engine_FU']+ftis['rh_engine_FU']
    t = ftis['time']
    W, M, X_cg = mcg(ftis['time'][(t==tstart)].item(), fuel_used[(t==tstart)], 1)
    ac = ac_B21(m=W/9.80665)

    # Phugoid
    print('PHUGOID')
    fig, ax = plt.subplots()
    Sym_2(tstart, tend, ftis, ac)

    syss = sym_flight(ac)
    T = ftis['time'][(t >= tstart) & (t<tend)]
    U = ftis['delta_e'][(t >= tstart) & (t<tend)]*np.pi/180
    # U = np.ones_like(U)#*np.average(U)
    X0 = np.array([0, 0, 0, 0])
    t, yout, xout = control.forced_response(syss, U=U, T=T, X0=X0)
    # yout[0, :] *= -4
    # yout[1, :] *= 3
    # yout[2, :] *= -2
    # yout[3, :] *= -3
    plot_response(T, yout[0, :], yout[1, :], yout[2, :], yout[3, :], "TAS [m/s]", "AOA [deg]", "Pitch angle [deg]", "Pitch rate [deg/s]")

    plt.tight_layout()
    plt.show()
    plt.plot()
    # plt.plot(T,U)
    # plt.show()

    # # Short Period
    # print('SHORT PERIOD')
    # Sym_2(52*60,53*60, ftis, ac)
    # plt.show()


    # # Dutch Roll
    # print('DUTCH ROLL')
    # tstart = 55 * 60  # [s]
    # tend = 57 * 60  # [s]
    # fuel_used = ftis['lh_engine_FU'] + ftis['rh_engine_FU']
    # t = ftis['time']
    #
    # W, M, X_cg = mcg(ftis['time'][(t == tstart)].item(), fuel_used[(t == tstart)], 1)
    # ac = FlightParams(m=W / 9.80665)
    #
    #
    # fig, ax = plt.subplots()
    # Asym_2(tstart, tend, ftis, ac)
    #
    # syss = sym_flight(ac)
    # T = ftis['time'][(t >= tstart) & (t < tend)]
    # U = ftis['delta_e'][(t >= tstart) & (t < tend)] * np.pi / 180
    # # U = np.ones_like(U)#*np.average(U)
    # X0 = np.array([0, 0, 0, 0])
    # t, yout, xout = control.forced_response(syss, U=U, T=T, X0=X0)
    # yout[0, :] *= -4
    # yout[1, :] *= 3
    # yout[2, :] *= -2
    # yout[3, :] *= -3
    # plot_response(T, yout[0, :], yout[1, :], yout[2, :], yout[3, :], "TAS [m/s]", "AOA [deg]", "Pitch angle [deg]",
    #               "Pitch rate [deg/s]")
    #
    # plt.tight_layout()
    # plt.show()
    # plt.plot()





    # # Dutch Roll YD
    # Asym_1(51*60 +13, 52*60, ftis)
    # # Aperiodic Roll
    # print('APERIODIC ROLL')
    # Asym_1(53*60 +30,54*60, ftis)
    # print('SPIRAL')
    # Asym_1(60*60,65*60, ftis)

