# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:55:25 2020

Author: Sander Orbons

Eigenvalue finding
"""
from scipy.signal import find_peaks  
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


def eig(f_list,t_list):
    peaks, _ = find_peaks(f_list)
    t = []
    f = []
    for i in range(len(peaks)):
        t.append(t_list[peaks[i]]-t_list[peaks[0]])
        f.append(f_list[peaks[i]])
    def e(t,a,b,c):
        return a*np.exp(-b*t) + c
        
    popt, pcov = curve_fit(e,t,f)
   
    
    A_half = (f[0]-popt[2])/2
    curve = []
    for r in range(len(t_list)):
        tc = t_list[r]-t_list[peaks[0]]
        curve.append(e(tc,popt[0],popt[1],popt[2]))
        
    
    for j in range(len(t_list)):
         if curve[j] <= A_half:
                       
              T_half = t_list[j] - t_list[peaks[0]]
              break
           
    P = t[1]-t[0]    
    XI = np.log(0.5)/T_half
    ETA = (np.pi*2)/P
    return T_half, P, XI, ETA

#________________________________________________________________________________#
    

def ploteig(f_list, t_list):
    peaks, _ = find_peaks(f_list)
    t = []
    f = []
    for i in range(len(peaks)):
        t.append(t_list[peaks[i]]-t_list[peaks[0]])
        f.append(f_list[peaks[i]])
    def e(t,a,b,c):
        return a*np.exp(-b*t) + c
        
    popt, pcov = curve_fit(e,t,f)
   
    
    #curve = e(t_list,popt[0],popt[1],popt[2])
    #print('Curve', curve)
    A_half = (f[0]-popt[2])/2
    curve = []
    for r in range(len(t_list)):
        tc = t_list[r]-t_list[peaks[0]]
        curve.append(e(tc,popt[0],popt[1],popt[2]))
        
    
    for j in range(len(t_list)):
         if curve[j] <= A_half:
                       
              T_half = t_list[j] - t_list[peaks[0]]
              break
           
    P = t[1]-t[0]    
    XI = np.log(0.5)/T_half
    ETA = (np.pi*2)/P
    plt.plot(curve,t_list)
    plt.show()
    return T_half, P, XI, ETA
    