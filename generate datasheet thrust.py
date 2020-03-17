#make dataset for thrust
import numpy as np

rawtrim = [[8177,161,5,0,3,0,424,470,708,4],
           [8432,149,6,-0.5,3,-17,420,464,728,2],
           [8540,140,7,-0.9,3,-22.5,415,461,742,2],
           [8672,131,9,-1.45,3,-41.5,413,456,772,2],
           [7835,171,5,0.25,3,31.5,430,476,792,4],
           [7475,182,4,0.6,3,54,435,483,810,6],
           [7485,191,4,0.85,2.95,88,443,491,826,6]]

rawcl = [[9003.3333, 248.666666667, 1.75, 744.5, 795.5, 418.3333333, 9.066666667],
        [8996.6667, 218.833333333, 2.5, 639.1666667, 680.166666667, 457,7],
        [9010,190.1666667,3.716666667,496.8333333,535.8333333,535.5000,5.16666667],
        [9010,158.8333333,5.800000,435.0000,487.8333333,561.1666667,3.533333333],
        [9036.67,131.5000,8.700000,384.5000,419.0000,596.0000,2.5],
        [8956.67,118.0000,10.83333333,391.1666667,439.,624.8333333,1.8]]


pound = 0.45359237
hour  = 60*60
foot  = 0.3048
knots = 0.514444
T_0   = 288.15-273
l     = 0.0065
gamma = 1.4
g_0   = 9.881
R     = 287.058
p_0   = 101325
rho_0 = 1.225

dat = []

for i in range(6):
    h_p = rawcl[i][0]*foot
    
    p   = p_0*(1+ l*h_p/T_0)**(-g_0/(l*R))
    Vc  = rawcl[i][1]*knots
    M   = np.sqrt(2/(gamma-1)*((1+p_0/p*((1+(gamma-1)*rho_0*Vc**2/(2*gamma*p_0))**((gamma-1)/gamma)-1))**((gamma-1)/gamma)-1))
    
    TAT = rawcl[i][6]
    T_ISA  = T_0 + l*h_p
    T_static = TAT/(1+((gamma-1)/2*M**2))
    dTemp  = T_static - T_ISA
    TATcor = TAT - dTemp
    dT_ISA = TATcor - T_ISA

    MFL = rawcl[i][3]/hour*pound
    MFR = rawcl[i][4]/hour*pound

    dat.append([h_p,M,dT_ISA,MFL,MFR])

trimdat = []
for i in range(7):
    h_p = rawtrim[i][0]*foot

    p   = p_0*(1+ l*h_p/T_0)**(-g_0/(l*R))
    Vc  = rawtrim[i][1]*knots
    M   = np.sqrt(2/(gamma-1)*((1+p_0/p*((1+(gamma-1)*rho_0*Vc**2/(2*gamma*p_0))**((gamma-1)/gamma)-1))**((gamma-1)/gamma)-1))

    TAT = rawtrim[i][9]
    T_ISA  = T_0 + l*h_p
    T_static = TAT/(1+((gamma-1)/2*M**2))
    dTemp  = T_static - T_ISA
    TATcor = TAT - dTemp
    dT_ISA = TATcor - T_ISA

    MFL = rawtrim[i][6]/hour*pound
    MFR = rawtrim[i][7]/hour*pound

    dat.append([h_p,M,dT_ISA,MFL,MFR])

np.savetxt('matlab.dat',dat, delimiter = ' ')
