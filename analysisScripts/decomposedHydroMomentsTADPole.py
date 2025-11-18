from scipy.special import eval_legendre

import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from matplotlib.colors import LogNorm

from scipy.interpolate import griddata
from scipy import integrate

from tqdm import tqdm
import os
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

###########################################################
### Set up argparse
###########################################################
parser = argparse.ArgumentParser(
  description='Plot flow around TADPole, and extract dipolar component. Please include path to data.')
parser.add_argument('dataPath', type=str, help="Path to the data")
args = parser.parse_args()
root_directory = args.dataPath
target_filename = "/swimmerflowfield.dat"

# If you have repeats (recommended), modify the next line appropriately to take it into account. Will also need to look at lines 227, and the plotting at line 197.
filepath = [str(root_directory)]

print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")

# First loop reads data and averages
for files in tqdm(filepath):

    cnt=0

    # Extracts simulation input parameters using the json file 
    InputFile=files+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)

    L=tion['domain']
    dt=int(tion.get('dt',10))
    t=int(tion['simSteps']/dt)
    swt=tion["swFlowOut"]/dt
    dipole=tion['dsSwim']
    sigma=tion['sigSwim']
    len_dip=sigma*(dipole+0.5)*0.5


    if len(L)==2:
        Lx=tion['domain'][0]
        Ly=tion['domain'][1]
        Lz=1

    if len(L)==3:
        Lx=tion['domain'][0]
        Ly=tion['domain'][1]
        Lz=tion['domain'][2]

    Size=Lx*Ly*Lz

    X=np.zeros(Lx)
    Y=np.zeros(Ly)
    Z=np.zeros(Lz)

    vRav=np.zeros((Lx,Ly,Lz,3))

    file=files+'/swimmerflowfield.dat'

    infile = open(file,"r")
    # print( '\tToss header ...' )
    for i in range(13):
        #toss header
        line = infile.readline()
        #print line

    # print( '\tRead data ...' )
    i=0
    j=0
    n=-1
    while infile:
        
        line = infile.readline()
        if not line:
            break
        else:
            split = line.split("\t")
        
        if len(split) == 7:
            time,Qx,Qy,Qz,Vx,Vy,Vz = split # to take into account old vs new format
        else:
            break # error handling

        j=int(i%(Size))
        i=i+1

        kx=int((j//Lz)//Ly)
        ky=int((j//Lz)%Ly)
        kz=int(j%Lz)
        
        Qx=float(Qx);Qy=float(Qy);Qz=float(Qz)
        if Vx!='' and Vy!='' and Vz!='':
            Vx=float(Vx);Vy=float(Vy);Vz=float(Vz)
        else:
            Vx=0.0;Vy=0.0;Vz=0.0
        

        dx=Qx-Lx/2+0.5
        dy=Qy-Ly/2+0.5
        dz=Qz-Lz/2+0.5
        cnt+=1

        X[kx]=dx
        Y[ky]=dy
        Z[kz]=dz

        vRav[kx,ky,kz,0]+=Vx/swt
        vRav[kx,ky,kz,1]+=Vy/swt
        vRav[kx,ky,kz,2]+=Vz/swt

    np.save(files+'/x',X)
    np.save(files+'/y',Y)
    np.save(files+'/z',Z)
    np.save(files+'/v',vRav/cnt)

# Second loop performs an azimuthal average (in 2-D, this is just a folding-over of the top and the bottom half)
for files in tqdm(filepath):

    Exp=files
    X=np.load(Exp+'/x.npy')
    Y=np.load(Exp+'/y.npy')
    Z=np.load(Exp+'/z.npy')

    V=np.load(Exp+'/v.npy')

    InputFile=Exp+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)

    L=tion['domain']

    if len(L)==3:
        diags=int(np.ceil(len(Y)/np.sqrt(2)))

    if len(L)==2:
        diags=int(np.ceil(len(Y)/2))

    xyU=np.zeros((len(X),diags,3))
    cnt=np.zeros((len(X),diags))

    ix=-1
    for x in X:
        ix+=1
        iy=-1
        for y in Y:
            iy+=1
            iz=-1
            for z in Z:
                iz+=1

                phi=0-np.arctan2(z,y)
                y2=y*np.cos(phi)-z*np.sin(phi)
                z2=y*np.sin(phi)+z*np.cos(phi)

                y2=int(y2+len(Y)/2)
                z2=int(z2+len(Z)/2)

                while y2>=len(Y):
                    y2-=len(Y)
                while y2<0:
                    y2+=len(Y)
                
                y2-=int(len(Y)/2)

                while y2<0:
                    y2+=len(Y)
            
                xyU[ix,y2,0]+=V[ix,iy,iz,0]
                xyU[ix,y2,1]+=V[ix,iy,iz,1]*np.cos(phi)-V[ix,iy,iz,2]*np.sin(phi)
                xyU[ix,y2,2]+=V[ix,iy,iz,1]*np.sin(phi)+V[ix,iy,iz,2]*np.cos(phi)
                cnt[ix,y2]+=1

    offset=-tion['sigSwim']*tion['dsSwim']
    offset=-tion['sigSwim']-tion['dsSwim']
    offset=int(tion['sigSwim']*(1-2*tion['dsSwim'])/4)

    a=xyU[:,:,0]/cnt
    b=xyU[:,:,1]/cnt
    c=xyU[:,:,2]/cnt

    for i in range(len(a[:,0])):
        j=i+offset
        if j<0:
            j+=L[0]
        xyU[j,:,0]=a[i,:]
        xyU[j,:,1]=b[i,:]
        xyU[j,:,2]=c[i,:]
    np.save(str(files)+'/xyU',xyU)

    # Plotting the flow field (do this after averaging over repeats, if you are going to do this)

    plt.figure(figsize=(10, 6), dpi=80)

    for i in range(len(xyU[:,0,0])):
            for j in range(int(len(xyU[0,:,0]))):
                if (i-int(L[0]/2))*(i-int(L[0]/2))+j*j>=L[0]*L[0]/4:
                    xyU[i,j,0]=0
                    xyU[i,j,1]=0

    xyU=xyU[:,0:int(L[1]/2),:]

    xyU2=np.sqrt(xyU[:,:,0]*xyU[:,:,0]+xyU[:,:,1]*xyU[:,:,1]).T


    u,v=np.meshgrid(np.arange(len(xyU[:,0,0])),np.arange(len(xyU[0,:,0])))
    plt.xlabel('Distance along swimmer axis, $r$',fontsize=32)
    plt.ylabel(r'Distance away\\ from swimmer axis',fontsize=32)

    skip=1

    w2=np.sqrt(xyU[:,:,0]*xyU[:,:,0]+xyU[:,:,1]*xyU[:,:,1])
    w=[min(i[0],np.power(10.0,-15)) for i in w2]

    a2=plt.imshow(xyU2,norm=LogNorm(),label='Vel')
    plt.quiver(u[::skip,::skip],v[::skip,::skip],xyU[::skip,::skip,0].T,-xyU[::skip,::skip,1].T,linewidth=w)
    plt.savefig(str(files)+'/FlowAroundTADPole.pdf')

# Extracting the dipolar component of the flow, then plotting its radial evolution.
for files in tqdm(filepath):
    xyU=np.load(files+'/xyU.npy')

    InputFile=files+'/../input.json'
    with open(InputFile, 'r') as InputF:
        tion=json.load(InputF)
    L=tion['domain']

    if len(L)==3:
        dim=3
        Lx=L[0]
        Ly=int(L[1]/np.sqrt(2))

    if len(L)==2:
        dim=2
        Lx=L[0]
        Ly=int(L[1]/2)

    # Construct Ur (different from xyU2, dot product)
    # Also creating data structures for the interpolation
    Ur=[]
    grid=[]

    for i in range(Lx):
        for j in range(Ly):
            x=float(i-L[0]/2)+0.5
            y=1*(float(j)+0.5)

            Ur.append((xyU[i,j,0]*x+xyU[i,j,1]*y)/np.sqrt(x*x+y*y))
            grid.append([x,y])
        
    for i in range(Lx):
        for j in range(Ly):
            x=float(i-L[0]/2)+0.5
            y=-1*(float(j)+0.5)

            Ur.append((xyU[i,j,0]*x-xyU[i,j,1]*y)/np.sqrt(x*x+y*y))
            grid.append([x,y])

    a=1
    n=3

    # Creating the array of radii that will be considered (staying close to the swimmer for decent stats)
    rl=[0.5+i*a for i in range(int(np.ceil(Lx/2)/a)-1)]

    # The data structure for each U_{r,n}(r)
    url=np.zeros((len(rl),n))

    er=-1
    d=Ly

    def Interp(x,y):
        return griddata(grid,Ur,(x,y),method='cubic')

    N=300

    dt=np.pi/N
    Th=np.zeros((N,n))

    for i in range(N):
        for k in range(n):
            if dim==2:
                Th[i,k]=np.cos(dt*k*i)
            if dim==3:
                Th[i,k]=np.sin(dt*i)*eval_legendre(k,np.cos(dt*i))

    # Loop over radii
    for r in rl:

        er+=1
        
        dt=np.pi/N
        th_list=np.zeros((N,n))
        tint=np.zeros(N)

        X=[]
        Y=[]

        # Loop over polar shells
        for ti in (range(N)):
            th=ti*dt
            
            tint[ti]=th
            # Finding the x,y position within the loop 
            c=np.cos(th)
            s=np.sin(th)
            x=r*c
            y=r*s
            X.append(x)
            Y.append(y)
            A=1

        for k in (range(n)):
            th_list[:,k]=griddata(grid,Ur,(X,Y),method='cubic')
            th_list[:,k]*=Th[:,k]
            url[er,k]=integrate.simpson(th_list[:,k],tint,dt)
            url[er,k]*=(0.5*(2.0*k+1.0))


    # Plotting the result

    # # Version with n moments 
    # for k in range(n):
    #     plt.plot(rl,abs(url[:,k]),label=k)

    # Dipolar moment only version 
    lim=int(len_dip)
    plt.cla()
    plt.figure(figsize=(10, 6), dpi=80)
    plt.plot(rl[0:lim+1],abs(url[0:lim+1,2]),linestyle='--',color=ed.teal)
    plt.plot(rl[lim:-1],abs(url[lim:-1,2]),linestyle='-',color=ed.teal)
    plt.scatter(rl[lim],abs(url[lim,2]),marker='o',color=ed.crimson,label='Theoretical size of dipole')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radial distance from TADPole centre')
    plt.ylabel(r'Intensity of dipolar \\component of TADPole flow')
    plt.savefig(files+'/DipolarFlow.pdf')

    np.save(files+'/url',url)

    for k in range(n):
        np.save(files+'/url'+str(k),url[:,k])
