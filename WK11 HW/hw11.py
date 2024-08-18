# -*- coding: utf-8 -*-
"""
ASTE-546 EM_PIC Example, Part B
"""

import numpy
import pylab as pl


#computes curl of E
def computeCurlE(Ex, Ey):
    curl = numpy.zeros((ni-1,nj-1))
    for i in range(ni-1):
        for j in range(nj-1):
            S = 0.5*(Ex[i][j]+Ex[i+1][j])
            E = 0.5*(Ey[i+1][j]+Ey[i+1][j+1])
            N = -0.5*(Ex[i+1][j+1]+Ex[i][j+1])
            W = -0.5*(Ey[i][j+1]+Ey[i][j])
            curl[i][j] = (S*dx + E*dy + N*dx + W*dy)/(dx*dy)
    return curl

#computes curl B on [1:nx-1],[1-ny-1]
def computeCurlB(Bz):
    curlX = numpy.zeros((ni,nj))
    curlY = numpy.zeros((ni,nj))
    for i in range(ni-2):
        for j in range (nj-2):
            N = 0.5*(Bz[i+1][j+1] + Bz[i][j+1])
            S = 0.5*(Bz[i][j] + Bz[i+1][j])
            E = 0.5*(Bz[i+1][j] + Bz[i+1][j+1])
            W = 0.5*(Bz[i][j+1] + Bz[i][j])
            curlX[i+1][j+1] = (N-S)/dy
            curlY[i+1][j+1] = -(E-W)/dx
    return curlX,curlY

#plots data with a given title
def plot(ax,data,title=""):
    pl.sca(ax)
    pl.cla()
    cf = pl.contourf(pos_x, pos_y, numpy.transpose(data),8,alpha=.75,cmap='jet')
    ax.set_yticks(pos_y)
    ax.set_xticks(pos_x)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    pl.xlim(min(pos_x),max(pos_x))
    pl.ylim(min(pos_y),max(pos_y))
    #ax.grid(b=True,which='both',color='#444',linestyle='-')
    ax.set_aspect('equal', adjustable='box')
    pl.title(title)
    pl.colorbar(cf,ax=pl.gca(),orientation='vertical',shrink=0.75, pad=0.01)

#sets field to zero along domain boundaries to prevent reflection
def clearBoundaries(F):
    ni,nj = numpy.shape(F)
    F[0,:] = 0
    F[ni-1,:] = 0
    F[:,0] = 0
    F[:,nj-1] = 0
    
#constants
eps0 = 8.854187e-12     #permitivity of free space
mu0 = 1.2566370e-6      #permeability of free space
c2 = 1/(eps0*mu0)   #speed of light squared
c = numpy.sqrt(c2)  #speed of light
E = 1.602e-19     #elementary charge

#define electric mesh
ni = 81
nj = 41
dx = 0.025
dy = 0.025

#time step
dt = 0.5*(dx+dy)/c

Ex = numpy.zeros((ni,nj))
Ey = numpy.zeros((ni,nj))
den_i = numpy.zeros((ni,nj))

#define magnetic mesh
Bz = numpy.zeros((ni-1,nj-1))

#placeholder for current density
den = numpy.zeros((ni,nj))
den[30:50,18:22] = 1e4
rho = den*E

ux = 1e3
Jx = rho*ux
Jy = numpy.zeros((ni,nj))

#main loop
for it in range(100):
    dB_dt = -computeCurlE(Ex,Ey)
    Bz += dB_dt*dt
    
    curlBx, curlBy = computeCurlB(Bz)
    dEx_dt = c2*curlBx - Jx/eps0
    dEy_dt = c2*curlBy - Jy/eps0
    Ex += dEx_dt*dt
    Ey += dEy_dt*dt
    
    #clear fields on boundaries to prevent reflection
    clearBoundaries(Ex)
    clearBoundaries(Ey)
    clearBoundaries(Bz)


    
    
#plotting
fig = pl.figure(num=None, figsize=(12, 10), dpi=80, facecolor='w', edgecolor='k')
sub = (pl.subplot(3,2,1),pl.subplot(3,2,2),pl.subplot(3,2,3),pl.subplot(3,2,4),
       pl.subplot(3,2,5),pl.subplot(3,2,6))

pos_x = numpy.linspace(0,(ni-1)*dx,ni)
pos_y = numpy.linspace(0,(nj-1)*dy,nj)
plot(sub[0],Ex,title="Ex")
plot(sub[1],Ey,title="Ey")
plot(sub[3],rho,title="rho")
plot(sub[4],Jx,title="Jx")
plot(sub[5],Jy,title="Jy")

#recompute pos_x/pos_y for the magnetic mesh
pos_x = numpy.linspace(0.5*dx,(ni-2)*dx,ni-1)
pos_y = numpy.linspace(0.5*dy,(nj-2)*dy,nj-1)
plot(sub[2],Bz,title="Bz")



