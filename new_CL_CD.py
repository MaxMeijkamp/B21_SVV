import numpy as np
import xlrd as xl
import matplotlib.pyplot as plt

from flight_conditions import *
from MCG import mcg
from functools import partial


def get_cla_cl0(alphacl):
    alpha = alphacl[:,0]
    cl = alphacl[:,1]
    poly_clalpha = np.polyfit(alpha, cl, 1)
    cla = poly_clalpha[0]
    cl0 = poly_clalpha[1]
    return cla, cl0


def get_CD0_e(clcd, A):
    cl = clcd[:,0]
    cl2 = cl*cl
    cd = clcd[:,1]
    poly_cl2_cd = np.polyfit(cl2, cd, 1)
    cd0 = poly_cl2_cd[1]
    pi_e_A = 1 / poly_cl2_cd[0]
    e = pi_e_A / np.pi / A
    return cd0, e


def resultarrays(data):
    mcg2 = np.vectorize(partial(mcg, situation=0, unit=1))
    W, M, xcg = mcg2(data[:,3])
    ac = FlightParams()
    alpha_CL = []
    CL_CD = []
    for m, meas in zip(W/9.80665, data):
        ac.hp, ac.V0, ac.alpha0, ac.T0, ac.m = meas[0], meas[1], meas[2], meas[4], m
        alpha_CL.append([meas[2], ac.CL])
        CD = meas[5]/(0.5*1.225*ac.V*ac.V*ac.S)
        CL_CD.append([ac.CL, CD])
    return np.asarray(alpha_CL), np.asarray(CL_CD), ac.A


def get_data():
    # format of data: each row is a measurement point. columns represent:
    # altitude [m], IAS [m/s], alpha [rad], F_used [lbs], Temp [K], T [N]
    stationarytest = np.zeros((6, 8))
    stationarytestsheet = xl.open_workbook('Post_Flight_Datasheet_Flight_B21.xlsx').sheet_by_index(0)
    for i in range(6):
        for j in range(7):
            stationarytest[i][j] = stationarytestsheet.cell_value(27 + i, 3 + j)
    stationarytest[:,0] *= 0.3048
    stationarytest[:,1] *= 0.514444444
    stationarytest[:,2] *= np.pi/180.
    stationarytest[:,6] += 273.15
    thrustperengine = np.array([[2362.34, 2649.15],
                                [2009.85, 2235.36],
                                [1409.47, 1629.51],
                                [1262.78, 1579.92],
                                [1139.45, 1354.89],
                                [1281.91, 1595.86]])  # Thrust in left and right engine for stationary measurements
    for i in range(len(thrustperengine)):  # total thrust during stationary measurements
        stationarytest[i][7] = sum(thrustperengine[i])
    stationarytest = np.delete(stationarytest, (3,4), 1)

    reftest = np.zeros((6, 8))
    refsheet = xl.open_workbook('Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx').sheet_by_index(0)
    for i in range(6):
        for j in range(7):
            reftest[i][j] = refsheet.cell_value(27 + i, 3 + j)
    reftest[:,0] *= 0.3048
    reftest[:,1] *= 0.514444444
    reftest[:,2] *= np.pi/180.
    reftest[:,6] += 273.15

    thrustrefperengine = np.array([[3700.82,	3807.6],
    [3013.99,	3076.02],
    [2410.22,	2536.92],
    [1865.59,	2018.17],
    [1892.67,	2076.08],
    [2194.11,	2389.75]])  # Thrust in left and right engine for stationary measurements
    for i in range(len(thrustrefperengine)):  # total thrust during stationary measurements
        reftest[i][7] = sum(thrustrefperengine[i])
    reftest = np.delete(reftest, (3,4), 1)
    return stationarytest, reftest

if __name__ == "__main__":
    stattest, reftest = get_data()
    print("REFERENCE DATA")
    alpha_CL, CL_CD, A = resultarrays(reftest)
    cd0, e = get_CD0_e(CL_CD, A)
    cla, cl0 = get_cla_cl0(alpha_CL)
    print("CD0 = {}; e = {}".format(cd0, e))
    print("CLa = {}; CL0 = {}".format(cla, cl0))
    print()
    print("OWN DATA")
    alpha_CL, CL_CD, A = resultarrays(stattest)
    cd0, e = get_CD0_e(CL_CD, A)
    cla, cl0 = get_cla_cl0(alpha_CL)
    print("CD0 = {}; e = {}".format(cd0, e))
    print("CLa = {}; CL0 = {}".format(cla, cl0))