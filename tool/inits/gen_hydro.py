from pyFC import LogNormalFractalCube, write_cube
import pyFC
#import matplotlib.pyplot as pl
#import matplotlib.cm as cm

Nx    = 256 
Ny    = 256 
Nz    = 256 
kmin  = 10.
mean  = 1.
sigma = 5**0.5
beta  = -5./3.

Lx    = 50.
Ly    = 50.
Lz    = 50.


fc = LogNormalFractalCube(ni=Nx, nj=Ny, nk=Nz, kmin=kmin, kmax=None, mean=mean, sigma=sigma, beta=beta)

fc.gen_cube()

write_cube(fc=fc, fname='Hydro_IC', app=True, prec='single',boxsize=[Lx,Ly,Lz])
