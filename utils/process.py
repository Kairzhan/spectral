# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>. 
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from numpy import ma
from matplotlib import ticker, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import cos, sin

import argparse
parser = argparse.ArgumentParser(description='Postprocessing of binary Fortran files.')
parser.add_argument('-v' , '--visc', help='Viscosity coefficient', required=False, type=float)
parser.add_argument('-n' , '--size', help='Mesh resolution (cubic)', required=False, type=int)
parser.add_argument('-t' , '--dt', help='Timestep', required=False, type=float)
parser.add_argument('-s' , '--start', help='Start step', required=False, type=int)
parser.add_argument('-e' , '--end', help='End step', required=False, type=int)
parser.add_argument('-i' , '--step', help='Iteration step', required=False, type=int)
parser.add_argument('-l' , '--levels', help='Number of contour levels', required=False, type=int)
parser.add_argument('-m' , '--mode', help='Mode. Compute contours, energy, skewness? Example (switch off all): -m 0 0 0', required=False, nargs='+', type=int)
# parser.add_argument('-c','--contours', help='Draw contours', required=False, action=argparse.BooleanOptionalAction)
# parser.add_argument('-k','--contours', help='Compute KE, DKE', required=False, action=argparse.BooleanOptionalAction)
# parser.add_argument('-w','--contours', help='Compute Skewness and Kurtosis', required=False, action=argparse.BooleanOptionalAction)
args=parser.parse_args()

# Start/end time steps
istart=50
iend=2000
step=istart

# Mesh resolution (cubic mesh is assumed)
n=64

# Viscosity and timestep
visc=0.003333432
dt=0.005

nlevels=30

# Plot contours? energy? skewness?
lcontours=True
lenergy=True
lskewness=True

# Substitute by cmd line parameters
pvisc=args.visc
if pvisc:
    visc=pvisc

pnx=args.size
if pnx:
    n=pnx

pdt=args.dt
if pdt:
    dt=pdt

pstart=args.start
if pstart:
    istart=pstart

pend=args.end
if pend:
    iend=pend

pstep=args.step
if pstep:
    step=pstep

plevels=args.levels
if plevels:
    nlevels=plevels

pmode=args.mode
if pmode:
    lcontours=bool(pmode[0])
    lenergy=bool(pmode[1])
    lskewness=bool(pmode[2])

# print(nlevels,lcontours, lenergy, lskewness)
# print(f"viscosity: {visc}")
# print(f"nx: {n}")
# print(f"nlevels: {nlevels}")
# print(f"start: {istart}")
# print(f"end: {iend}")
# print(f"step: {step}")
# print(f"mode: {lcontours} {lenergy} {lskewness}")
# print(f"dt: {dt}")
# #print(f": {}")

# Velocity field
ux=np.zeros(n*n*n)
uy=np.zeros(n*n*n)
uz=np.zeros(n*n*n)

# Wavenumbers
kx=np.zeros(n*n*n).reshape(n,n,n)
ky=np.zeros(n*n*n).reshape(n,n,n)
kz=np.zeros(n*n*n).reshape(n,n,n)

for j in range(n):
    for k in range(n):
        kx[:,j,k]=np.fft.fftfreq(n, 1./n)

for i in range(n):
    for k in range(n):
        ky[i,:,k]=np.fft.fftfreq(n, 1./n)

for i in range(n):
    for j in range(n):
        kz[i,j,:]=np.fft.fftfreq(n, 1./n)

k2=kx**2+ky**2+kz**2
kabs=np.sqrt(k2)

# Half wavenumbers (for r2c transforms)
kxh=np.fft.fftfreq(n, 1/n)
kyh=np.fft.fftfreq(n, 1/n)
kzh=np.fft.fftfreq(n, 1/n)
KX, KY, KZ=np.meshgrid(kxh[:], kyh[:], kzh[:n//2+1], indexing='ij')

# Remove previous postprocessing
if lenergy and os.path.exists("stat_te_spectral.dat"):
    os.remove("stat_te_spectral.dat")
    
if lenergy and os.path.exists("stat_te.dat"):
    os.remove("stat_te.dat")
    
if lskewness and os.path.exists("skew-flat.dat"):
    os.remove("skew-flat.dat")

# Time loop
for istep in range(istart, iend+step, step):
    print(istep)
    
    ux=np.fromfile(f'uuu{istep:09}.dat', dtype=np.float32)
    uy=np.fromfile(f'vvv{istep:09}.dat', dtype=np.float32)
    uz=np.fromfile(f'www{istep:09}.dat', dtype=np.float32)

    ux=ux.reshape(n,n,n, order='F')
    uy=uy.reshape(n,n,n, order='F')
    uz=uz.reshape(n,n,n, order='F')
    
    # Contours
    if lcontours:
        x=np.linspace(0, 2*np.pi, n)
        y=np.linspace(0, 2*np.pi, n)
        X, Y = np.meshgrid(x, y)

        # Setup contour levels first
        if istep==istart:
            levelsx=np.linspace(np.amin(ux), np.amax(ux), nlevels)
            levelsy=np.linspace(np.amin(ux), np.amax(uy), nlevels)
            #levelsz=np.linspace(np.amin(ux), np.amax(uz), nlevels)
            levelsz=np.linspace(-1, 1, nlevels)
            print(levelsx, levelsy, levelsz)

        # UX
        Z=ux[n//2,:,:].reshape(n,n).tolist()
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)
        cs = ax.contourf(X, Y, Z, levels=levelsx,cmap=cm.coolwarm) #PuBu_r)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cs, cax=cax)
        plt.savefig(f'Contours/Contour-U-{istep:09}.png', dpi=120)
        plt.close()

        # UY
        Z=uy[n//2,:,:].reshape(n,n).tolist()
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)
        cs = ax.contourf(X, Y, Z, levels=levelsy, cmap=cm.coolwarm) #PuBu_r)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cs, cax=cax)
        plt.savefig(f'Contours/Contour-V-{istep:09}.png', dpi=120)
        plt.close()

        # UZ
        Z=uz[n//2,:,:].reshape(n,n).tolist()
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)
        cs = ax.contourf(X, Y, Z, levels=levelsz, cmap=cm.coolwarm) #PuBu_r)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cs, cax=cax)
        plt.savefig(f'Contours/Contour-W-{istep:09}.png', dpi=120)
        plt.close()
        
    # Kinetic energy
    if lenergy:
        ke=0.5*np.sum(ux*ux+uy*uy+uz*uz)/n**3

        # Kinetic energy from DFT
        ux_hat=np.fft.fftn(ux)/n**3
        uy_hat=np.fft.fftn(uy)/n**3
        uz_hat=np.fft.fftn(uz)/n**3

        et=0.5*np.real(ux_hat*np.conj(ux_hat)+uy_hat*np.conj(uy_hat)+uz_hat*np.conj(uz_hat))
        ets=2*visc*et*k2
        
        ke_spectral=np.sum(np.where(kabs<n/2, et, 0))
        ed_spectral=np.sum(np.where(kabs<n/2, ets, 0))
        
        with open("stat_te.dat", "a") as f:
            f.write("{} {} {}\n".format(istep, istep*dt, ke))
    
        with open("stat_te_spectral.dat", "a") as f:
            f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(istep, istep*dt, ke_spectral, ed_spectral, 0, 0, 0, 0, 0, 0, 0))

        # Compute spectrum
        ett=np.zeros(n)
        et=np.fft.fftshift(et)
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    wn=int(np.rint(np.sqrt((i-n/2)**2+(j-n/2)**2+(k-n/2)**2+0.0)))
                    ett[wn]=ett[wn]+et[i,j,k]

        with open(f'spe-{istep:09}.dat', "w") as f:
            for wn in range(n):
                f.write("{} {}\n".format(wn, ett[wn]))
