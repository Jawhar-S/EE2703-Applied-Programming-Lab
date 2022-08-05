'''
------------------------------------------------------------------------------------------------------
 EE2703 - Assignment-5
 Done by Jawhar S (EE20B049)
 Last modified on: 07/03/2022
 Description: This program solves for the potential at the surface of the copper plate by solving the 
    Laplace Equation for a given set of parameters. This is done by employing discrete differentiation. 
    The current density at the surface of the copper plate is also obtained.
------------------------------------------------------------------------------------------------------
'''
from pylab import *
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from argparse import ArgumentParser 

# Function to find ideal coefficients A, B for the fit y = A*exp(Bx)
def lstsq_fit(x,y):

    x = x.reshape((-1,1))
    one_arr = np.ones((x.shape[0],1))
    xvar = concatenate((one_arr,x), axis=1)
    logA, B = lstsq(xvar, log(y), rcond=None)[0]
    return exp(logA), B

# Parsing the arguments if given any
arg_def = ArgumentParser()

arg_def.add_argument('--Nx', default=25, type=int, help='X Resolution')
arg_def.add_argument('--Ny', default=25, type=int, help='Y Resolution')
arg_def.add_argument('--radius', default=8, type=int,help='Radius in indices')
arg_def.add_argument('--Niter', default=1500, type=int, help='Number of Iterations')

argv = arg_def.parse_args()
Nx, Ny, Nr, Niter = argv.Nx, argv.Ny, argv.radius, argv.Niter

# Initializing potential array
phi = zeros((Ny,Nx), dtype='float64')

x = linspace(0,Nx-1,Nx) - Nx//2             
y = linspace(0,Ny-1,Ny) - Ny//2             
Y, X = meshgrid(y,x)
ii = np.where(X*X + Y*Y <= Nr*Nr)
phi[ii] = 1.0                              
ii_T = np.array(ii).T

# Contour plot with 1V points marked red
figure(1)
contourf(Y,X,phi)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')            
title("Contour plot of Potential")
xlabel(r"x $\rightarrow$")
ylabel(r"y $\rightarrow$")
show()

oldphi = phi.copy()
errors = np.zeros(Niter)
# Performing iteration
for k in range(Niter):
    oldphi = phi.copy() 
    # Updating potential array
    phi[1:-1, 1:-1] = 0.25*(oldphi[1:-1, 0:-2] + oldphi[1:-1, 2:] + oldphi[0:-2, 1:-1] + oldphi[2:, 1:-1])
    # Boundary conditions
    phi[1:-1, 0] = phi[1:-1,1]          # Left 
    phi[1:-1, -1] = phi[1:-1, -2]       # Right 
    phi[0, 1:-1] = phi[1, 1:-1]         # Top 
    phi[-1, 1:-1] = 0                   # Bottom (Ground)
    # Re-Updating the wire potential
    phi[ii] = 1.0                       
    # Obtaining the maximum error vector 
    errors[k] = (abs(phi-oldphi)).max()

# Fitting the error
ind = array(range(0,Niter))
ind1 = array(range(0,Niter))
ind2 = ind1[500:]
# fit with all points
A1,B1 = lstsq_fit(ind1, errors)    
# fit with all except first 500
A2,B2 = lstsq_fit(ind2, errors[500:])        
fit1 = A1*exp(B1*ind1)
fit2 = A2*exp(B2*ind2)

# semilog plots
figure(2)
semilogy(ind[::50],errors[::50], 'ro--', markersize=12)
semilogy(ind1[::50], fit1[::50], 'go--', markersize=8)
semilogy(ind2[::50], fit2[::50], 'yo--', markersize=4)
legend(('error', 'fit all', 'fit 500+'))
title("Error and Error in Fitted Plots")
xlabel(r"Iteration Count $\rightarrow$")
ylabel(r"Error $\rightarrow$")
show()

# loglog plot
figure(3)
loglog(ind[::50],errors[::50], 'ro--', markersize=8)
title("Error Plot (LogLog)")
xlabel(r"Iteration Count $\rightarrow$")
ylabel(r"Error $\rightarrow$")
show()

# Error estimation
maxerr_1 = abs(-A1*exp(B1*(Niter+0.5)))
maxerr_2 = abs(-A2*exp(B2*(Niter+0.5)))

print("Error Bound for fit 1:", maxerr_1)
print("Error Bound for fit 2:", maxerr_2)

# Surface Plot
fig4 = figure(4)
ax = p3.Axes3D(fig4)
ax.plot_surface(Y,X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
xlabel(r"Y $\rightarrow$")
ylabel(r"X $\rightarrow$")
ax.set_zlabel(r"Potential $\rightarrow$")
ax.set_title("Surface Plot of Potential")
show()

# Contour Plot
figure(5)
contourf(Y,X[::-1], phi)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')                
title("Contour Plot of Potential")
xlabel(r"X $\rightarrow$")
ylabel(r"Y $\rightarrow$")
show()

# Obtaining Current density at the surface 
Jx = np.zeros(phi.shape)
Jy = np.zeros(phi.shape)
Jx[:,1:-1] = (phi[:,0:-2]-phi[:,2:])/2
Jy[1:-1,:] = (phi[2:, :]-phi[0:-2,:])/2

# Vector Plot of Current 
figure(6)
quiver(Y,X[::-1],Jx,Jy, scale=6)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')    
title("Vector Plot of Current Density")
xlabel(r"X $\rightarrow$")
ylabel(r"Y $\rightarrow$")
show()
