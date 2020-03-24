# Author: Merlijn Hunik
# C_L vs C_D calculation

import math as m
import numpy as np
import xlrd as xl

from flight_conditions import *
from MCG import mcg
import matplotlib.pyplot as plt


import scipy.interpolate as sc

R = 287.058
g = 9.80665
S = 30
Ws = 14815.95*g*0.453592
A = 30 / 2.0569 ** 2


# input a matrix with [h[ft], Vias[kt], a[deg], FFl[lbs/hr], FFr[lbs/hr], , F.used[lbs], TAT[C], T[N]]
# if plot = 0 a cdCl^2 plot will be produced, otherwise a plot of cl alpha will be made

def cdclalfa(stationarytest, plot):
    R = 287.058
    g = 9.80665
    S = 30
    Ws = 14815.95 * g * 0.453592
    rho0 = 1.225

    for i in range(len(stationarytest)):
        stationarytest[i][8] = ISA_calculator(stationarytest[i][0] * 0.3044, R, g)[1]

    C_l = np.zeros([6,1])
    C_D = np.zeros([6,1])
    alfa= np.zeros([6,1])
    C_lsquare = np.zeros([6,1])

    for i in range(len(stationarytest)):
        W = mcg(stationarytest[i][5], 0, 1)[0]
        rho = stationarytest[i][8]
        #V0 = stationarytest[i][1] * 0.51444 * np.sqrt(rho/rho0)*np.sqrt(Ws/W)
        V0 = stationarytest[i][1]* 0.514444
        T = stationarytest[i][7]
        hp = stationarytest[i][0]*0.3044
        C_l[i][0] = FlightParams(m=W/g, V0=V0, hp=hp,  rho0=rho, T0=(stationarytest[i][7])+273.15).CL
        #C_l[i][0] = 2*W/rho/V0/V0/S
        C_D[i][0] = C_l[i][0]*T/W
        alfa[i][0] = stationarytest[i][2]/180*np.pi
        C_lsquare[i][0] = C_l[i][0]**2

    corrcl2cd = np.polyfit(C_lsquare.T[0], C_D.T[0], 1)
    trendlinecdcl = np.poly1d(corrcl2cd)

    corrC_l_alfa = np.polyfit(alfa.T[0], C_l.T[0], 1)
    trend_alfa_cl = np.poly1d(corrC_l_alfa)

    C_l_alfa = corrC_l_alfa[0]
    Cl0 = corrC_l_alfa[1]
    print('Cl_alpha = ', C_l_alfa, ', Cl0=', Cl0)

    Cd0 = corrcl2cd[1]
    pi_e_A = 1 / corrcl2cd[0]
    e = pi_e_A / np.pi / A

    print('Cdo =', Cd0,", e =", e)

    if plot == 0:
        plt.plot(C_lsquare, C_D, 'o', C_lsquare, trendlinecdcl(C_lsquare))
        plt.title(r'$C_D$-$C_L^2$ curve')
        plt.xlabel(r"$C_L^2$ [-]")
        plt.ylabel("$C_D$ [-]")
    else:
        plt.plot(alfa, C_l, 'o', alfa, trend_alfa_cl(alfa))
        plt.title(r'$C_l$-$\alpha$ curve')
        plt.xlabel(r"$\alpha$ [rad]")
        plt.ylabel("$C_l$ [-]")

    plt.show()
    return C_l, C_D, alfa, C_lsquare

stationarytest = np.zeros((6, 9))
stationaryexcell = xl.open_workbook('Post_Flight_Datasheet_Flight_B21.xlsx')
stationarytestsheet = stationaryexcell.sheet_by_index(0)
# for loop to give an array with 6x [h[ft], V[kt], a[deg], FFl[lbs/hr], FFr[lbs/hr] F.used[lbs], TAT[C], T[N], rho[kg/m^3]
for i in range(6):
    for j in range(7):
        stationarytest[i][j] = stationarytestsheet.cell_value(27 + i, 3 + j)

thrustperengine = np.array([[2362.34, 2649.15],
                            [2009.85, 2235.36],
                            [1409.47, 1629.51],
                            [1262.78, 1579.92],
                            [1139.45, 1354.89],
                            [1281.91, 1595.86]])  # Thrust in left and right engine for stationary measurements
for i in range(len(thrustperengine)):  # total thrust during stationary measurements
    stationarytest[i][7] = sum(thrustperengine[i])

reftest = np.zeros((6, 9))
refexcell = xl.open_workbook('Post_Flight_Datasheet_Flight_1_DD_12_3_2018.xlsx')
refsheet = refexcell.sheet_by_index(0)
# for loop to give an array with 6x [h[ft], V[kt], a[deg], FFl[lbs/hr], FFr[lbs/hr] F.used[lbs], TAT[C], T[N], rho[kg/m^3]
for i in range(6):
    for j in range(7):
        reftest[i][j] = refsheet.cell_value(27 + i, 3 + j)

thrustrefperengine = np.array([[3700.82,	3807.6],
[3013.99,	3076.02],
[2410.22,	2536.92],
[1865.59,	2018.17],
[1892.67,	2076.08],
[2194.11,	2389.75]])  # Thrust in left and right engine for stationary measurements
for i in range(len(thrustrefperengine)):  # total thrust during stationary measurements
    reftest[i][7] = sum(thrustrefperengine[i])

a = cdclalfa(reftest, 1)






