#!/usr/bin/env python

import numpy as np
import healpy as hp
import random

def genRMmapQU(alpha, nside,seed, sigma):
    ell=np.arange(3*nside,dtype=float)
    cl=np.zeros(len(ell))
    beta=0.
    lmin=2
    for l in range(len(ell)):
        if l<=lmin: cl[l]=0
        else: cl[l]=10**beta*ell[l]**alpha
            
        #Cl_Q,U created
    count=0
    lmax=len(cl)
    alm_Q_wo=np.zeros(lmax*(lmax+1)/2,dtype=complex)
    alm_U_wo=np.zeros(lmax*(lmax+1)/2,dtype=complex)

    ##set seed 
    np.random.seed(seed)
    for l in range(lmax):
        alm_Q_wo[count:count+l+1]=np.sqrt(cl[l])/np.sqrt(2)*np.random.normal(0,1,l+1)+1j*np.sqrt(cl[l])/np.sqrt(2)*np.random.normal(0,1,l+1)
        alm_U_wo[count:count+l+1]=np.sqrt(cl[l])/np.sqrt(2)*np.random.normal(0,1,l+1)+1j*np.sqrt(cl[l])/np.sqrt(2)*np.random.normal(0,1,l+1)
        count+=l+1
    
    #rearrange alm
    alm_Q=np.zeros(lmax*(lmax+1)/2,dtype=complex)
    alm_U=np.zeros(lmax*(lmax+1)/2,dtype=complex)
    count=0
    for m in range(lmax):
        for l in range(m,lmax):
            alm_Q[count]=alm_Q_wo[l*(l+1)/2+m]
            alm_U[count]=alm_U_wo[l*(l+1)/2+m]
            count+=1
    map_Q=hp.sphtfunc.alm2map(alm_Q, nside)
    map_U=hp.sphtfunc.alm2map(alm_U, nside)
    map_QU=np.zeros((len(map_Q),2))
    ###note maps are resclaled
    map_QU[:,0]=map_Q/np.std(map_Q)*sigma
    map_QU[:,1]=map_U/np.std(map_U)*sigma
    
    
    return map_QU
