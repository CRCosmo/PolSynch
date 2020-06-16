#!/usr/bin/env python

import numpy as np
import healpy as hp
import random

#function to generate the NxNxlmax matrix to use for correlations
def genRho(nRM,lmax,rho_info, corr):
    rho=np.zeros((nRM,nRM,lmax))
    if corr==True: 
        print "Reading correlation info"
        for rmi in range(nRM):
            for rmj in range(nRM):
                for l in range(lmax):
                    if l<200: rho[rmi,rmj,l]=(rho_info[0]*l+rho_info[1])
                    else: rho[rmi,rmj,l]=(rho_info[2]*l**rho_info[3])
    return rho

#function to create covariance matrix to use in RM map generation. If rho is zero no correlation are added.
def RM_cov(nRM, ell, alpha, rho):
    #nRM number N of RM values
    #alpha is a N vector
    #rho is a NxNxlmax matrix
    cov=np.zeros((nRM,nRM),dtype=float)
    for i in range(nRM):
        for j in range(i,nRM):
            if (j==i): cov[i,i]=ell**alpha[i]
            else: 
                cov[i,j]=0.75*np.sqrt(rho[i,j,ell])*np.sqrt(ell**alpha[i]*ell**alpha[j]) #rho is under squared root for Q,U to obtain rho in P
                cov[j,i]=cov[i,j]
    return cov  

def genRMcube_corr(nRM, nside,seed, cov,sigma):
    #sigma is N vector
    #alpha is N vector
    #cov is (lmax,nRM,nRM)
    
    ell=np.arange(3*nside,dtype=float)

        #Cl_Q,U created
    count=0
    lmax=3*nside
    alm_Q_wo=np.zeros((lmax*(lmax+1)/2,nRM),dtype=complex)
    alm_U_wo=np.zeros((lmax*(lmax+1)/2,nRM),dtype=complex)

    ##set seed 
    np.random.seed(seed)
    for l in range(lmax):
            if l<=2: alm_Q_wo[count:count+l+1,:]=alm_U_wo[count:count+l+1,:]=0
            else:
                alm_Q_wo[count:count+l+1,:]=1./np.sqrt(2)*(np.random.multivariate_normal(np.zeros(nRM),cov[l,:],l+1)+1j*np.random.multivariate_normal(np.zeros(nRM),cov[l,:],l+1))
                alm_U_wo[count:count+l+1,:]=1./np.sqrt(2)*(np.random.multivariate_normal(np.zeros(nRM),cov[l,:],l+1)+1j*np.random.multivariate_normal(np.zeros(nRM),cov[l,:],l+1))
            count+=l+1
    
    #rearrange alm
    alm_Q=np.zeros((lmax*(lmax+1)/2,nRM),dtype=complex)
    alm_U=np.zeros((lmax*(lmax+1)/2,nRM),dtype=complex)
    count=0
    for m in range(lmax):
        for l in range(m,lmax):
            alm_Q[count,:]=alm_Q_wo[l*(l+1)/2+m,:]
            alm_U[count,:]=alm_U_wo[l*(l+1)/2+m,:]
            count+=1
    cubeQ=np.zeros((nRM,hp.nside2npix(nside)),dtype=float)
    cubeU=np.zeros((nRM,hp.nside2npix(nside)),dtype=float)
    cube_QU=np.zeros((nRM,hp.nside2npix(nside),2))
    for rm in range(nRM):
        alm_Q_C=alm_Q[:,rm]
        alm_U_C=alm_U[:,rm]
        cubeQ[rm,:]=hp.sphtfunc.alm2map(np.ascontiguousarray(alm_Q_C), nside)
        cubeU[rm,:]=hp.sphtfunc.alm2map(np.ascontiguousarray(alm_U_C), nside)
        
        ###note maps are resclaled
        cube_QU[rm,:,0]=cubeQ[rm,:]/np.std(cubeQ[rm,:])*sigma[rm]
        cube_QU[rm,:,1]=cubeU[rm,:]/np.std(cubeU[rm,:])*sigma[rm]
    
    
    return cube_QU
