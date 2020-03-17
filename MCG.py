# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:10:12 2020

@author: Sander Orbons

Mass and x_c.g. determination
"""
import numpy as np


def mcg(w_pax, m_f, unit):
    # w_pax: np array of weight passengers and crew, seat 1->10          [pounds]
    # m_f: fuel weight in pounds at given time                       [pounds]
    # unit: 0-> output in imperial units, 1-> output in metric units 
    
    # Empty weight 
    m_BEM = 9165.0
    x_BEM = 291.65
    # Seat locations
    x_seat = np.array([131,131,214,214,251,251,288,288,170])
    # Fuel c.g. arm
    list_x_f = np.array([298.16, 591.18, 879.08, 1165.42, 1448.40, 1732.53, 2014.80, 2298.84, 2581.92, 2866.30, 3150.18, 3434.52, 3718.52, 4003.23, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.90, 5994.04, 6278.47, 6562.82, 6846.96, 7131.00, 7415.33, 7699.60, 7984.34, 8269.06, 8554.05, 8839.04, 9124.80, 9410.62, 9696.97, 9983.40, 10270.08, 10556.84, 10843.87, 11131.00, 11418.20, 11705.50, 11993.31, 12281.18, 12,569.04, 12856.86, 13144.73, 13432.48, 13720.56, 14008.46, 14320.34])
    x_f = list_x_f[int(m_f/100)]
    # compute Moment and Weight
    M = sum(w_pax*x_seat+m_BEM*x_BEM+m_f*x_f)
    W = sum(w_pax)+m_BEM+m_f 
    # convert to metric units if needed [kg & m]
    if unit == 1:
        M = M*0.0254
        W = W*0.453592
    # compute x_c.g.
    X_cg = M/W
    print('mass and center of gravity have been computed')
    return W, M, X_cg