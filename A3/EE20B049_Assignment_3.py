'''
----------------------------------------------------------------------------------------------------------------
 EE2703 - Assignment-3
 Done by Jawhar S (EE20B049)
 Last modified on: 18/02/2022
 Description: This program extracts data from 'fitting.dat'; plots the data along with the noise levels and the exact
  function; calculates the MS error over discrete data points and identifies the minima; plots MS error vs
  different standard deviation of noises on both linear and loglog scales.
----------------------------------------------------------------------------------------------------------------
'''
import numpy as np
import scipy.special as scp
from pylab import *

def g(t,A,B):
    g=A*scp.jn(2,t)+B*t
    return g

# Extracting the data from 'fitting.dat' and storing them in arrays 
dataCols=[]
dataCols=np.loadtxt('fitting.dat',dtype=float)
t=np.array(dataCols[:,0])            
f=np.asarray(dataCols)[:,1:]  
# Array consisting of values uniformly sampled on a logarthmic scale
STDS=np.logspace(-1,-3,9)                         
# Plotting the data along with the true curve
figure(figsize=(7,6))  
for i in range(9):
    plot(t,f[:,i],label=r'$\sigma$=%.3f'%STDS[i])
plot(t,g(t,1.05,-0.105),'-k',label='True curve')
title("Q4: Data to be fitted in theory",fontsize=16)   
legend()                             
xlabel(r't  $\rightarrow$',fontsize=14)                 
ylabel(r'f(t) + noise  $\rightarrow$',fontsize=14)        
grid(True) 
show()                                           

# Plotting an Errorbar of every 5th data item along with the true curve 
data=f[:,0]
figure(figsize=(7,6)) 
errorbar(t[::5],data[::5],STDS[0],fmt='ro',label='Errorbar') 
plot(t,g(t,1.05,-0.105),'k',label='f(t)') 
title("Q5: Data points for $\sigma$ = 0.10 along with exact function",fontsize=16)                 
legend()                                     
xlabel(r't  $\rightarrow$',fontsize=14)                         
grid(True)  
show()                                                  

# Asserting Equality of g(t,A0,B0) and M.p
x=scp.jn(2,t)
M=c_[x,t]
p=np.array([1.05,-0.105])
LHS=np.matmul(M,p)      
RHS=np.array(g(t,1.05,-0.105))
assert(np.allclose(LHS,RHS))    

# Plotting the Contour and the Minima
err_matrix=np.zeros([21,21])           
A=np.linspace(0,2,21)               
B=np.linspace(-0.2,0,21)            
col1=f[:,0]
for i in range(21):
    for j in range(21):
        # Mean Square Error
        err_matrix[i][j] = np.sum((col1-np.array(g(t,A[i],B[j])))**2)/101  
figure(figsize=(7,6)) 
Contour=contour(A,B,err_matrix,20) 
title(r"Q8: Contour plot of $\epsilon_{ij}$",fontsize=16)
xlabel(r'A  $\rightarrow$',fontsize=14)                 
ylabel(r'B  $\rightarrow$',fontsize=14)                 
clabel(Contour,Contour.levels[0:5],inline=1,fontsize=10)  
# Obtaining the location of Minima     
Min=np.unravel_index(np.argmin(err_matrix),err_matrix.shape)
plot(A[Min[0]],B[Min[1]],'or',markersize=5)                               
annotate('(%0.2f,%0.2f)'%(A[Min[0]],B[Min[1]]),(A[Min[0]],B[Min[1]]))     
grid(False)                   
show() 

# Plotting MSE vs Standard deviation on a linear scale
figure(figsize=(7,6)) 
MSE=[np.linalg.lstsq(M,f[:,i],rcond=None)[0] for i in range(9)]                     
MSE=np.asarray(MSE)                                                        
A_Err = abs(MSE[:,0]-1.05)                                                 
B_Err = abs(MSE[:,1]+0.105)                                               
plot(STDS,A_Err,'ro--',label='Aerr')                                
plot(STDS,B_Err,'go--',label='Berr')                                
legend()
title("Q10: Variation of error with noise",fontsize=16)                                              
xlabel(r'Noise standard deviation  $\rightarrow$',fontsize=14)          
ylabel(r'MS error$\rightarrow$',fontsize=14)                          
grid(True)  
show()                                                       

# Plotting MSE vs Standard deviation on a loglog scale
figure(figsize=(7,6))                                 
loglog(STDS,A_Err,'ro',label='Aerr')                
stem(STDS,A_Err,'-ro', use_line_collection=True)                              
loglog(STDS,B_Err,'go',label='Berr')                
stem(STDS,B_Err,'-go',use_line_collection=True)                            
legend()     
title("Q11: Variation of error with noise",fontsize=16)                      
xlabel(r'$\sigma_{n}  \rightarrow$',fontsize=14)        
ylabel(r'MS error$  \rightarrow$',fontsize=14)          
grid(False)                                          
show() 





