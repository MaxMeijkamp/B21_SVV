# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:10:12 2020
@author: Sander Orbons
Mass and x_c.g. determination
"""
import numpy as np


def mcg(m_f, situation, unit):
    # m_f:       fuel used in pounds at given time [pounds]
    # situation: 1-> no passengers, 0-> everyone in correct seat, else-> one passenger to the front
    # unit:      0-> output in imperial units [lbs,lbsinch and inch], else-> output in metric units (N,Nm and m) 
    
    # Start weight 
    m_BEM = 9165.0
    x_BEM = 291.65
    #m_bf = 4100
    m_bf = 4050
    
    # cabine configuration
    #x_seat = np.array([131,131,170,150,214,214,251,251,288,288])
    x_seat = np.array([131,131,170,134,214,214,251,251,288,288])
    if situation == 0:
       # m_pax = 2.20462*np.array([80,102,68,0,76,70,73,82,62,90.5])
        m_pax = 2.20462*np.array([95,92,74,0,66,61,75,78,86,68])
    elif situation == 1:
        m_pax = np.zeros(10)
    else:
        #m_pax = 2.20462*np.array([80,102,68,90.5,76,70,73,82,62,0])
        m_pax = 2.20462*np.array([95,92,74,86,66,61,75,78,0,68])

    # Fuel Moment as from Appendix E
    list_M_f = np.array([298.16, 591.18, 879.08, 1165.42, 1448.40, 1732.53, 2014.80, 2298.84, 2581.92, 2866.30, 3150.18, 3434.52, 3718.52, 4003.23, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.90, 5994.04, 6278.47, 6562.82, 6846.96, 7131.00, 7415.33, 7699.60, 7984.34, 8269.06, 8554.05, 8839.04, 9124.80, 9410.62, 9696.97, 9983.40, 10270.08, 10556.84, 10843.87, 11131.00, 11418.20, 11705.50, 11993.31, 12281.18, 12,569.04, 12856.86, 13144.73, 13432.48, 13720.56, 14008.46, 14320.34])
    M_f = list_M_f[int((-m_f+m_bf)/100)] 
    
    # compute Moment and Weight
    # W = lbs +
    W = sum(m_pax)+m_BEM+m_bf-m_f 
    M = sum(m_pax*x_seat) + m_BEM*x_BEM+M_f*100
    
    # convert to metric units if needed [N & m]
    if unit != 0:
        W = W*4.44822           #[lbs to N]
        M = M*0.11298482933333  #[lbsinch to Nm]
    
    # compute x_c.g.
    X_cg = M/W
    return W, M, X_cg