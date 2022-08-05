'''
----------------------------------------------------------------------------------------------------------------
 EE2703 - Assignment-4
 Done by Jawhar S (EE20B049)
 Last modified on: 26/02/2022
 Description: This program fits two functions, exp(x) and cos (cos(x)) over the interval [0, 2π) using:
    (i) exact formula of the coefficients in the fourier series 
    (ii) estimated coefficients of the fourier series by least squares method
----------------------------------------------------------------------------------------------------------------
'''

import numpy as np
from pylab import *
from scipy.integrate import quad

def exp(x):
    return np.exp(x)

# Periodic extension of the exponential function in the range [0, 2π)
def periodic_exp(x):
    return np.exp(x%(2*np.pi))

def coscos(x):
    return np.cos(np.cos(x))

# Cosine term that will be used in the calculation of fourier coefficients
def u(x, k, f):
    return f(x)*np.cos(k*x)

# Sine term that will be used in the calculation of fourier coefficients
def v(x, k, f):
    return f(x)*np.sin(k*x)

func_map={'exp(x)':exp,'cos(cos(x))':coscos}

# Function to calculate first 'n' fourier coefficients
def get_coefficients(n,func):
    coefficients=np.zeros(n)                                         
    func_str=func_map[func]
    # DC coefficient
    coefficients[0]=quad(func_str,0,2*np.pi)[0]/(2*np.pi) 
    # Cosine coefficients        
    for k in range(1,n,2):
        coefficients[k]=quad(u,0,2*np.pi,args=((k+1)/2, func_str))[0]/np.pi   
    # Sine coefficients
    for k in range(2,n,2):
        coefficients[k]=quad(v,0,2*np.pi,args=(k/2, func_str))[0]/np.pi        
    return coefficients

x=np.linspace(-2*np.pi,4*np.pi,400)        

# Semi-Log plot of exp(x) and periodic extension of exp(x) in the range [0, 2π)
figure(figsize=(7,6))                                         
semilogy(x,exp(x),'-r',label='Actual')            
semilogy(x,periodic_exp(x),'-g',label='Periodic Extension') 
title(r'$exp(x)$ and periodic extension of $exp(x)$',fontsize=16)     
xlabel(r'x$\rightarrow$',fontsize=14)                         
ylabel(r'$e^{x}\rightarrow$',fontsize=14)                     
legend(loc='upper right')                                     
grid(True)                                                    
show()

# Plot of cos(cos(x))
figure(figsize=(7,6))                                        
plot(x,coscos(x),'-r')    
title(r'cos(cos(x))',fontsize=16) 
xlabel(r'x$\rightarrow$',fontsize=14)                         
ylabel(r'$\cos(\cos(x))\rightarrow$',fontsize=14)             
grid(True)                                                    
show()

exp_coeffs=get_coefficients(51,'exp(x)')               
coscos_coeffs=get_coefficients(51,'cos(cos(x))')        

# Log-Log plot of fourier coefficients of exp(x)
figure(figsize=(7,6))                                        
loglog(range(51),np.abs(exp_coeffs),'ro')                
xlabel(r'n$\rightarrow$',fontsize=14)                         
ylabel(r'Magnitude of the coefficients$\rightarrow$',fontsize=14)     
title(r'Fourier Coefficients of $e^{x}$',fontsize=16)         
grid(True)                                                    
show()

# Semi-Log plot of fourier coefficients of exp(x)
figure(figsize=(7,6))                                        
semilogy(range(51),np.abs(exp_coeffs),'ro')              
xlabel(r'n$\rightarrow$',fontsize=14)                         
ylabel(r'Magnitude of the coefficients$\rightarrow$',fontsize=14)     
title(r'Fourier Coefficients of $e^{x}$',fontsize=16)         
grid(True)                                                    
show()

# Log-Log plot of fourier coefficients of cos(cos(x))
figure(figsize=(7,6))                                        
loglog(range(51),np.abs(coscos_coeffs),'ro')             
xlabel(r'n$\rightarrow$',fontsize=14)                         
ylabel(r'Magnitude of the coefficients$\rightarrow$',fontsize=14)     
title(r'Fourier Coefficients of $cos(cos(x))$',fontsize=16)   
grid(True)                                                    
show()

# Semi-Log plot of fourier coefficients of cos(cos(x))
figure(figsize=(7,6))                                        
semilogy(range(51),np.abs(coscos_coeffs),'ro')           
xlabel(r'n$\rightarrow$',fontsize=14)                         
ylabel(r'Magnitude of the coefficients$\rightarrow$',fontsize=14)     
title(r'Fourier Coefficients of $cos(cos(x))$',fontsize=16)   
grid(True)                                                    
show()

# Least Squares method of estimating the fourier coefficients
x=np.linspace(0,2*np.pi,401)
x=x[:-1]
y=np.linspace(0,2*np.pi,400)
A=np.zeros((400,51))
A[:,0]=1
for k in range(1,26):
    A[:,2*k-1]=np.cos(k*x)
    A[:,2*k]=np.sin(k*x) 
b_exp=exp(x)  
b_coscos=coscos(x)
c_exp=np.linalg.lstsq(A,b_exp,rcond=None)[0]
c_coscos=np.linalg.lstsq(A,b_coscos,rcond=None)[0]

# Log-Log plot of fourier coefficients of exp(x) estimated by Least Square method
figure(figsize=(7,6))   
loglog(range(51),np.abs(c_exp),'go',label='Estimated Value')          
loglog(range(51),np.abs(exp_coeffs),'ro',label='True Value')     
xlabel(r'n$\rightarrow$',fontsize=14)                                 
ylabel(r'$Coefficient\rightarrow$',fontsize=14)                       
legend(loc='upper right')                                             
title(r'Fourier Coefficients Of $e^{x}$',fontsize=16)                 
grid(True)                                                            
show()

# Semi-Log plot of fourier coefficients of exp(x) estimated by Least Square method
figure(figsize=(7,6))  
semilogy(range(51),np.abs(c_exp),'go',label='Estimated Value')         
semilogy(range(51),np.abs(exp_coeffs),'ro',label='True Value')    
xlabel(r'n$\rightarrow$',fontsize=14)                                 
ylabel(r'$Coefficient\rightarrow$',fontsize=14)                       
legend(loc='upper right')                                             
title(r'Fourier Coefficients Of $e^{x}$',fontsize=16)                 
grid(True)                                                            
show()

# Log-Log plot of fourier coefficients of cos(cos(x)) estimated by Least Square method
figure(figsize=(7,6))  
loglog(range(51),np.abs(c_coscos),'go',label='Estimated Value')       
loglog(range(51),np.abs(coscos_coeffs),'ro',label='True Value')  
xlabel(r'n$\rightarrow$',fontsize=14)                                 
ylabel(r'$Coefficient\rightarrow$',fontsize=14)                       
legend(loc='upper right')                                             
title(r'Fourier Coefficients Of $cos(cos(x))$',fontsize=16)           
grid(True)                                                            
show()

# Semi-Log plot of fourier coefficients of cos(cos(x)) estimated by Least Square method
figure(figsize=(7,6))  
semilogy(range(51),np.abs(c_coscos),'go',label='Estimated Value')      
semilogy(range(51),np.abs(coscos_coeffs),'ro',label='True Value') 
xlabel(r'n$\rightarrow$',fontsize=14)                                 
ylabel(r'$Coefficient\rightarrow$',fontsize=14)                       
legend(loc='upper right')                                             
title(r'Fourier Coefficients Of $cos(cos(x))$',fontsize=16)           
grid(True)                                                            
show()

# Maximum Deviation between the fourier coefficients estimated by Least Squares method
max_exp_deviation=np.max(np.abs(exp_coeffs-c_exp))
max_coscos_deviation=np.max(np.abs(coscos_coeffs-c_coscos))

print("Largest Deviation in fourier coefficients of exp(x):",max_exp_deviation)
print("Largest Deviation in fourier coefficients of cos(cos(x)):",max_coscos_deviation)

# Plotting exp(x) using the estimated fourier coffieicents by Least Squares method 
exp_estimated=np.dot(A,c_exp)                                     
figure(figsize=(7,6)) 
semilogy(x,exp_estimated,'go',label='Estimated value')           
semilogy(x,exp(x),'-r',label='True value')
title(r'$exp(x)$ using the estimated fourier coffieicents',fontsize=16)                  
xlabel(r'n$\rightarrow$',fontsize=14)                     
ylabel(r'$e^{x}\rightarrow$',fontsize=14)                 
legend(loc='upper right')                                 
grid(True)                                                
show()

# Plotting cos(cos(x)) using the estimated fourier coffieicents by Least Squares method 
coscos_estimated=np.dot(A,c_coscos)                              
figure(figsize=(7,6))  
plot(x,coscos_estimated,'go',label='Estimated value')            
plot(x,coscos(x),'-r',label='True value')  
title(r'$cos(cos(x))$ using the estimated fourier coffieicents',fontsize=16)                 
xlabel(r'n$\rightarrow$',fontsize=14)                     
ylabel(r'$\cos(\cos(x))\rightarrow$',fontsize=14)         
legend(loc='upper left')                                  
grid(True)                                                
show()