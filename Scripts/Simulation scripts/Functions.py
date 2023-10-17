# -*- coding: utf-8 -*-
"""
Functions to define interaction parameter values
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands

"""
#%reset -f
#%%
import os
import sys

#some functions to define interaction parameters
def set_z(y,n,g,r,f):
    n=5
    g=-3/60
    f=1.65
    d = 0.2/60# natural death rate k_d0
#    r = maxgrowthrate # SPECIFIC FOR EACH SPECIES!!! and with n also specifc per drug
    if g != None: # for cidal drugs
        E_C = 1+ (((r-g)/(d)-1)*(f)**n)/((f)**n + ((g)/(d-r)))
        E_C_y = 1+ (((r-g)/(d)-1)*(f/(1+y))**n)/((f/(y+1))**n + ((g)/(d-r)))
        z = (d*(E_C_y-E_C)-r)/(-r) - 1
    elif g == None: # for static drugs
        E_S = 1 - ((f)**n)/(((f)**n) + d/(r-d))
        E_S_y = 1 - (((f)/(y+1))**n)/((((f)/(y+1))**n) + d/(r-d)) 
        z = (E_S_y)/(E_S) - 1
    return z

def set_x(y):
    x=-y/(y+1)
    return x

def scale_intradius(v1,r1,r2): #taking value v1 for radius r1 and calculating value v2 for r2
    N_C=6.8
    v2 = (2*r1+1)**2 / (2*r2+1)**2 *v1
    return v2
    
#Strength of interactions with a small radius (r=0) & & $f_{r2}=\frac{(2r_{1}+1)^{2}}{(2r_{1}+1)^{2}}f_{r1}$  & The strength of interactions with a small radius (r=0), are determined from interactions with a large radius (r=1). This formula is used to adjust the strength of interactions for different radii. The strength of the interaction for radius r2 is determined from the strength of the interaction at radius r1 multiplied with the ratio of available patches in each scenario. To convert interaction strength from r=1 to r=0, $f_{r=0}=9*f_{r=1}$
