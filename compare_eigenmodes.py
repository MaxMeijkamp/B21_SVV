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
from Cit_par_ref import *
from Cit_par_final import *


def get_ftis_data(filename, ref=False):
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

    if ref:
        keys = np.delete(keys, 9)

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


def Asym_2(t_start, t_end, ftis, ac: FlightParams):

    # Gathering data
    t = ftis['time']
    TAS_list = ftis['Dadc1_tas'][(t>=t_start) & (t<t_end)]
    roll_list = ftis['Ahrs1_Roll'][(t>=t_start) & (t<t_end)]*np.pi/180
    rollr_list = ftis['Ahrs1_bRollRate'][(t>=t_start) & (t<t_end)]*np.pi/180
    yawr_list = ftis['Ahrs1_bYawRate'][(t>=t_start) & (t<t_end)]*np.pi/180
    roll_list = roll_list - np.average(roll_list)
    rollr_list = (rollr_list-np.average(rollr_list))*ac.b*0.5/TAS_list[0]
    yawr_list = (yawr_list-np.average(yawr_list))*ac.b*0.5/TAS_list[0]
    t_list = t[(t>=t_start) & (t<t_end)]

    # Plotting figures for angle response
    plot_response(t_list, TAS_list, roll_list, rollr_list, yawr_list, "TAS [m/s]", "Roll angle [rad]", "Roll rate [rad/s]", "Yaw rate [rad/s]")

#________________________________________________________________#
    
# Symmetric response

def Sym_2(t_start,t_end, ftis, ac: FlightParams):

    t = ftis['time']
    aoa_list = ftis['vane_AOA'][(t>=t_start) & (t<t_end)]
    pitch_list = ftis['Ahrs1_Pitch'][(t>=t_start) & (t<t_end)]
    pitchr_list = ftis['Ahrs1_bPitchRate'][(t >= t_start) & (t < t_end)]
    TAS_list = ftis['Dadc1_tas'][(t>=t_start) & (t<t_end)]
    pitchr_list = np.deg2rad(pitchr_list) * ac.c / TAS_list[0] # Reduce and make dimensionless
    aoa_list = np.deg2rad(aoa_list - np.average(aoa_list)) # Reduce and change units
    pitch_list = np.deg2rad(pitch_list - np.average(pitch_list)) # Reduce and change units
    # pitch_list = np.deg2rad(pitch_list - pitchr_list[0]) # Reduce and change units
    TAS_list = (TAS_list-np.average(TAS_list))/TAS_list[0] # Make dimensionless and reduce
    t_list = t[(t>=t_start) & (t<t_end)]

    # Plotting figures for angle response

    plot_response(t_list, TAS_list, aoa_list, pitch_list, pitchr_list, "TAS [m/s]", "AOA [rad]", "Pitch angle [rad]", "Pitch rate [rad/s]")
    return


if __name__ == "__main__":
    # # ftis, units = get_ftis_data('FTISxprt-20200311_flight3.mat')
    # ftis, units = get_ftis_data('matlab.mat', ref=True)
    # tstart = 53*60+58 # [s]
    # tend = 57*60 # [s]
    # fuel_used = ftis['lh_engine_FU']+ftis['rh_engine_FU']
    # t = ftis['time']
    # W, M, X_cg = mcg(fuel_used[(t==tstart)].item(), 0,  1)
    # ac = ac_fin(m=W/9.80665)
    #
    # # Phugoid
    # print('PHUGOID')
    # fig, ax = plt.subplots()
    # Sym_2(tstart, tend, ftis, ac)
    #
    # syss = sym_flight(ac)
    #
    # T = ftis['time'][(t >= tstart) & (t<tend)]
    # U = ftis['delta_e'][(t >= tstart) & (t<tend)]*np.pi/180
    # # U = np.ones_like(U)#*np.average(U)
    # # X0 = np.array([0, 0, ftis['Ahrs1_Pitch'][(t==tstart)].item()*np.pi/180,
    # #                ftis['Ahrs1_bPitchRate'][(t==tstart)].item()*ac.c/ftis['Dadc1_tas'][(t==tstart)].item()*np.pi/180])
    #
    # X0 = np.array([0, 0, 0, ftis['Ahrs1_bPitchRate'][(t == tstart)].item() * ac.c / ftis['Dadc1_tas'][
    #                    (t == tstart)].item() * np.pi / 180])
    #
    # # X0 = np.array([0, 0, 0, 0])
    # t, yout, xout = control.forced_response(syss, U=U, T=T, X0=X0)
    # # yout[0, :] *= -4
    # # yout[1, :] *= 3
    # # yout[2, :] *= -2
    # # yout[3, :] *= -3
    # plot_response(T, yout[0, :]-np.average(yout[0, :]), yout[1, :]-np.average(yout[1, :]), yout[2, :]-np.average(yout[2, :]), yout[3, :], "TAS [m/s]", "AOA [rad]", "Pitch angle [rad]", "Pitch rate [rad/s]")
    #
    # plt.tight_layout()
    # plt.show()

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
    # ftis, units = get_ftis_data('FTISxprt-20200311_flight3.mat')
    ftis, units = get_ftis_data('matlab.mat', ref=True)
    tstart = 61*60+57 # [s]
    tend = 62*60+20 # [s]
    fuel_used = ftis['lh_engine_FU']+ftis['rh_engine_FU']
    t = ftis['time']
    W, M, X_cg = mcg(fuel_used[(t==tstart)].item(), 0,  1)
    ac = ac_fin(m=W/9.80665)

    # Phugoid
    print('DUTCH ROLL')
    fig, ax = plt.subplots()
    Asym_2(tstart, tend, ftis, ac)

    sysa = asym_flight(ac)

    T = ftis['time'][(t >= tstart) & (t<tend)]
    U = np.array([ftis['delta_a'][(t >= tstart) & (t<tend)]*np.pi/180, ftis['delta_r'][(t >= tstart) & (t<tend)]*np.pi/180])
    # U = np.ones_like(U)#*np.average(U)
    # X0 = np.array([0, 0, ftis['Ahrs1_Pitch'][(t==tstart)].item()*np.pi/180,
    #                ftis['Ahrs1_bPitchRate'][(t==tstart)].item()*ac.c/ftis['Dadc1_tas'][(t==tstart)].item()*np.pi/180])

    X0 = np.array([0, 0, 0, 0])

    # X0 = np.array([0, 0, 0, 0])
    t, yout, xout = control.forced_response(sysa, U=U, T=T, X0=X0)
    # yout[0, :] *= -4
    # yout[1, :] *= 3
    # yout[2, :] *= -2
    # yout[3, :] *= -3

    plot_response(T, yout[0, :]-np.average(yout[0, :]), yout[1, :]-np.average(yout[1, :]), yout[2, :]-np.average(yout[2, :]), yout[3, :], "TAS [m/s]", "Roll angle [rad]", "Roll rate [rad/s]", "Yaw rate [rad/s]")

    plt.tight_layout()
    plt.show()
