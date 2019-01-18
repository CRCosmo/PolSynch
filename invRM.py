#!/usr/bin/env python

import numpy as np
import healpy as hp
import os,sys
import random
from astropy.io import fits
import scipy.optimize 
from scipy import stats

c = 299792458.0 #speed of light

def freq2wl2(frequencies):
    """
    Convert *frequencies* (in MHz) to wavelength squared (in
    m^2)
    """
    return (c/(frequencies*1e6))**2

def comp_phases(wavelength_squared, phi):
    """
    Compute the phase factor exp(-2j*phi*wavelength_squared).
    """
    return np.exp(2j*phi*wavelength_squared)


def P_l2(npix, q_rm, u_rm, phi_array, freq):

    """
    Convert Q,U RM maps into Frequency space (RM^-1) 
    """

    wl2       = freq2wl2(freq)
    p_map    = np.zeros(npix, dtype = complex)

    
    num      = len(phi_array)

    for phi in range(num):
            print('processing frame '+str(phi+1)+'/'+str(num))
            p_rm = q_rm[phi] + 1.0j*u_rm[phi]
            phase   = comp_phases(wl2, phi)
            p_map += p_rm*phase
           
    p_map=p_map/num
        
    return p_map


    
    

