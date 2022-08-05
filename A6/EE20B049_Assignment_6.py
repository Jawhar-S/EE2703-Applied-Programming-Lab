'''
------------------------------------------------------------------------------------------
 EE2703 - Assignment-6
 Done by Jawhar S (EE20B049)
 Last modified on: 27/03/2022
 Description: This assignment analyses Linear Time-invariant Systems with numerical tools 
    in Python. Given set of problems in the assignment sheet are solved in this program.
------------------------------------------------------------------------------------------
'''
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp 

# P1: Obtaining time response for a decay of 0.5
freq = 1.5
decay = 0.5
num = np.poly1d([1,decay])
denom = np.polymul([1,2*decay,decay**2 + freq**2],[1,0,2.25])
H = sp.lti(num,denom)
t,x = sp.impulse(H,None,np.linspace(0,50,501))

plt.figure(1)
plt.plot(t,x)
plt.title(r"P1: Time Response (for decay = 0.5)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t) $\rightarrow$")
plt.show()

# P2: Obtaining time response for a decay of 0.05
freq = 1.5
decay = 0.05
num = np.poly1d([1,decay])
denom = np.polymul([1,2*decay,decay**2 + freq**2],[1,0,2.25])
H = sp.lti(num,denom)
t,x = sp.impulse(H,None,np.linspace(0,50,501))

plt.figure(2)
plt.plot(t,x)
plt.title(r"P2: Time Response (for decay = 0.05)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t) $\rightarrow$")
plt.show()

# P3: Obtaining time response using transfer function for different frequencies
H = sp.lti([1],[1,0,2.25])
t = np.linspace(0,50,501)

i = 3
for freq in np.linspace(1.4,1.6,5):
    f = np.cos(freq*t)*np.exp(-decay*t)
    t,x,svec = sp.lsim(H,f,t) 
    plt.figure(i)
    plt.plot(t,x)
    plt.title(r"P3: Time Response (freq = %.3f rad/s, decay = 0.05)" % freq)
    plt.xlabel(r"t $\rightarrow$")
    plt.ylabel(r"x(t) $\rightarrow$")
    i+=1
    plt.show()

# P4: Solving for time responses for the given coupled spring
num = [1,0,2]
denom = [1,0,3,0]
H = sp.lti(num,denom)
t = np.linspace(0,20,201)
t,x = sp.impulse(H, None, t)

num = [2]
denom = [1,0,3,0]
H = sp.lti(num,denom)
t = np.linspace(0,20,201)
t,y = sp.impulse(H, None, t)

plt.figure(8)
plt.plot(t,x, 'r')
plt.plot(t,y, 'g')
plt.title(r"P4 : Solution For The Coupled Spring")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t),y(t) $\rightarrow$")
plt.legend(['x(t)','y(t)'])
plt.show()

# Q5: RLC Circuit Transfer Funiion
H = sp.lti([1],[1e-12, 1e-4, 1])
w,S,phi = H.bode()

# Magnitude Response
plt.figure(9)
plt.subplot(2,1,1)
plt.semilogx(w,S)
plt.title("P5: Bode Plot of the RLC Filter")
plt.xlabel(r"w $\rightarrow$")
plt.ylabel(r"Magnitude(H) (dB) $\rightarrow$")

# Phase Response
plt.subplot(2,1,2)
plt.semilogx(w,phi)
plt.xlabel(r"w $\rightarrow$")
plt.ylabel(r"Phase(H) (degrees) $\rightarrow$")
plt.show()

# Q6: Obtaining output for a given input in the RLC filter
n = int(1e5)
t = np.linspace(0,1e-2, n+1)
f = np.cos(1e3 * t) - np.cos(1e6 * t)
t,v_o,svec = sp.lsim(H,f,t)

# Milliseconds scale
plt.figure(10)
plt.plot(t,v_o)
plt.title(r"P6: Output Response of the RLC Filter on long timescale (10ms)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"$v_{o}(t)$ $\rightarrow$")
plt.show()

# Microseconds scale
plt.figure(11)
plt.plot(t[:int(30*n/(1e4))],v_o[:int(30*n/(1e4))])
plt.title(r"P6: Output Response of the RLC Filter on short timescale (30$\mu$s)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"$v_{o}(t)$ $\rightarrow$")
plt.show()