'''
------------------------------------------------------------------------------------------
 EE2703 - Assignment-7
 Done by Jawhar S (EE20B049)
 Last modified on: 08/04/2022
 Description: This assignment solves Linear Time-invariant systems using Sympy and Scipy.
  The magnitude response, step response and outputs for decaying sinusiods and mixed sinusiods
  were obtained for two different systems. (i.e Lowpass and Highpass)
------------------------------------------------------------------------------------------
'''
import sympy as sm 
import matplotlib.pyplot as plt 
import numpy as np
import scipy.signal as sp 

PI = np.pi

# Function to solve Low Pass Circuit in Laplace domain
def lowpass(R1,R2,C1,C2,G,Vi):
    A = sm.Matrix([
                    [0,0,1,-1/G] \
                    , [-1/(1+s*R2*C2),1,0,0] \
                    , [0,-G,G,1] \
                    , [-(1/R1)-(1/R2)-s*C1,1/R2,0,s*C1]
                 ])
    
    b = sm.Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A,b,V)

# Function to solve High Pass Circuit in Laplace domain
def highpass(R1,R3,C1,C2,G,Vi):
    A = sm.Matrix([
                    [0,0,1,-1/G] \
                    , [-(s*C2*R3)/(1+s*C2*R3),1,0,0] \
                    , [0,-G,G,1] \
                    , [-(1/R1)-s*C2-s*C1,s*C2,0,1/R1]
                 ])
    
    b = sm.Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)

# Converts Sympy expression to Scipy.signal compatible LTI system defintion
def sp_lti(expr):
    
    expr=sm.simplify(expr) 
    num,denom=sm.fraction(expr)                                                                
    num,denom=sm.Poly(num,s),sm.Poly(denom,s)                                 
    num,denom=num.all_coeffs(),denom.all_coeffs()                                   
    num,denom=[float(f) for f in num],[float(f) for f in denom]
    return sp.lti(num,denom)
    
# Q1: Low Pass Step Response   
s = sm.symbols('s')
A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1/s)
H = sp_lti(V[3])
t = np.linspace(0,1.5e-4,int(1e5)+1)
t, lp_step = sp.impulse(H,None,t)

plt.figure(1)
plt.plot(t,lp_step)
plt.title("Q1: Low Pass Step Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.show()

# Q2: High Pass Response to Mixed Sinusoids 
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
H = sp_lti(V[3])
t = np.linspace(0,0.00001,int(1e5)+1)
vi = np.sin(2e3 * PI * t) + np.cos(2e6 * PI * t)
t,hp_out,svec = sp.lsim(H,vi,t)

plt.figure(2)
plt.plot(t,hp_out)
plt.title("Q2: High Pass Response to Mixed Sinusoids")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.show()

# Q3: High Pass Magnitude Response
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
w = np.logspace(0,8,801)
ss = 1j*w
hf = sm.lambdify(s,Vo,'numpy')
v = hf(ss)

plt.figure(3)
plt.loglog(w,abs(v),lw=2)
plt.title("Q3: High Pass Magnitude Response")
plt.grid(True)
plt.xlabel(r'$\omega \rightarrow$')
plt.ylabel(r'$|H(j\omega)| \rightarrow$')
plt.show()

# Q4: High Pass Response to Damped Sinusoid
H = sp_lti(Vo)
t_low = np.linspace(0,1,int(1e5)+1)
t_high = np.linspace(0,1e-5,int(1e5)+1)
vi_high = np.sin(1e7 * PI * t_high)*np.exp(-1.5e5*t_high)
vi_low = np.sin(1e2 * PI * t_low)*np.exp(-1.5*t_low)
t_high,hp_hf_out,svec = sp.lsim(H,vi_high,t_high)
t_low,hp_lf_out,svec = sp.lsim(H,vi_low,t_low)

plt.figure(4)
plt.plot(t_low,hp_lf_out)
plt.title("Q4: Highpass Decaying Low-Freq-Sinusoidal Input Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.show()

plt.figure(5)
plt.plot(t_high,hp_hf_out)
plt.title("Q4: Highpass Decaying High-Freq-Sinusoidal Input Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.show()

# Q5: High Pass Step Response
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vo = V[3]
H_step = sp_lti(Vo)
t = np.linspace(0,0.001,int(1e5)+1)
t, hp_step = sp.impulse(H_step,None,t)

plt.figure(6)
plt.plot(t,hp_step)
plt.title("Q5: High Pass Step Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.show()
