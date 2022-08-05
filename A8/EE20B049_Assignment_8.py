'''
------------------------------------------------------------------------------------------
 EE2703 - Assignment-8
 Done by Jawhar S (EE20B049)
 Last modified on: 12/05/2022
 Description: This assignment is focussed on obtaining the DFT (Digital Fourier Transform)
 of periodic signals and the Gauss function
------------------------------------------------------------------------------------------
'''
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt

# Spectrum of sin(5t)
x = np.linspace(-2*np.pi,2*np.pi,257)
x=x[:-1]
y = np.sin(5*x)
Y = fft.fftshift(fft.fft(y))/256
w = np.linspace(-64,64,257)
w=w[:-1]
fig,fig1 = plt.subplots(2,1,figsize=[6,8])
fig1[0].plot(w,abs(Y))
fig1[0].set_ylabel("Magnitude")
fig1[0].set_xlim([-10,10])
fig1[1].plot(w,np.angle(Y),'ro')
fig1[1].set_ylabel("Phase")
fig1[1].set_xlim([-10,10])
ii = np.where(abs(Y)>1e-3)
fig1[1].plot(w[ii],np.angle(Y[ii]),'go')
fig.suptitle("Spectrum of sin(5t)")
fig1[0].grid("true")
fig1[1].grid("true")
plt.show()

# Spectrum of (1+0.1cos(t))cos(10t)
t = np.linspace(-4*np.pi,4*np.pi,513)
t = t[:-1]
y = (1 + 0.1*np.cos(t))*np.cos(10*t)
Y = fft.fftshift(fft.fft(y))/512
w = np.linspace(-64,64,513)
w = w[:-1]
fig,fig2 = plt.subplots(2,1,figsize=[6,8])
fig2[0].plot(w,abs(Y))
fig2[0].set_ylabel("Magnitude")
fig2[0].set_xlim([-15,15])
fig2[1].plot(w,np.angle(Y),'ro')
fig2[1].set_ylabel("Phase")
fig2[1].set_xlim([-15,15])
ii = np.where(abs(Y) > 1e-3)
fig2[1].plot(w[ii],np.angle(Y[ii]),'go')
fig.suptitle("Spectrum of (1+0.1cos(t))cos(10t)")
fig2[0].grid("true")
fig2[1].grid("true")
plt.show()


# Spectrum of sin^3(t)
x = np.linspace(-4*np.pi,4*np.pi,513)
x=x[:-1]
y = pow(np.sin(x),3)
Y = fft.fftshift(fft.fft(y))/512
w = np.linspace(-64,64,513)
w=w[:-1]
fig,fig3 = plt.subplots(2,1,figsize=[6,8])
fig3[0].plot(w,abs(Y))
fig3[0].set_ylabel("Magnitude")
fig3[0].set_xlim([-5,5])
fig3[1].plot(w,np.angle(Y),'ro')
fig3[1].set_ylabel("Phase")
fig3[1].set_xlim([-5,5])
ii = np.where(abs(Y)>1e-3)
fig3[1].plot(w[ii],np.angle(Y[ii]),'go')
fig.suptitle("Spectrum of sin(t)^3")
fig3[0].grid("true")
fig3[1].grid("true")
plt.show()

# Spectrum of cos^3(t)
t = np.linspace(-4*np.pi,4*np.pi,513)
t=t[:-1]
y = pow(np.cos(t),3)
Y = fft.fftshift(fft.fft(y))/512
w = np.linspace(-64,64,513)
w=w[:-1]
fig,fig4 = plt.subplots(2,1,figsize=[6,8])
fig4[0].plot(w,abs(Y))
fig4[0].set_ylabel("Magnitude")
fig4[0].set_xlim([-5,5])
fig4[1].plot(w,np.angle(Y),'ro')
fig4[1].set_ylabel("Phase")
fig4[1].set_xlim([-5,5])
ii = np.where(abs(Y)>1e-3)
fig4[1].plot(w[ii],np.angle(Y[ii]),'go')
fig.suptitle("Spectrum of cos(t)^3")
fig4[0].grid("true")
fig4[1].grid("true")
plt.show()

# Spectrum of cos(20*t+5*cos(t))
x = np.linspace(-4*np.pi,4*np.pi,513)
x=x[:-1]
y = np.cos(20*x + 5*np.cos(x))
Y = fft.fftshift(fft.fft(y))/512
w = np.linspace(-64,63,513)
w=w[:-1]
fig,fig5 = plt.subplots(2,1,figsize=[6,8])
fig5[0].plot(w,abs(Y))
fig5[0].set_ylabel("Magnitude")
fig5[0].set_xlim([-40,40])
fig5[1].set_ylabel("Phase")
fig5[1].set_xlim([-40,40])
ii = np.where(abs(Y)>1e-3)
fig5[1].plot(w[ii],np.angle(Y[ii]),'go')
fig.suptitle("Spectrum of cos(20t + 5 cos(t))")
fig5[0].grid("true")
fig5[1].grid("true")
plt.grid("true")
plt.show()

# Spectrum of Gauss Function (N=512)
t = np.linspace(-4*np.pi,4*np.pi,513)
t=t[:-1]
y = np.exp(-0.5*pow(t,2))
w = np.linspace(-64,64,513)
w=w[:-1]

Y = fft.fftshift(abs(fft.fft(y)))/512
# Normalise the gaussian transfer function
Y = Y*np.sqrt(2*np.pi)/max(Y)
Y2 = np.exp(-w**2/2)*np.sqrt(2*np.pi)
# Calculate error in estimate
err = abs(Y-Y2).max()

fig,fig6 = plt.subplots(2,1,figsize=[6,8])
fig6[0].plot(w,abs(Y2),w,abs(Y))
fig6[0].set_ylabel("Magnitude")
fig6[0].set_xlim([-10,10])
fig6[1].set_ylabel("Phase")
fig6[1].set_xlim([-10,10])
ii = np.where(abs(Y)>1e-2)
fig6[1].plot(w[ii],np.unwrap(np.angle(Y[ii])),'go')
fig.suptitle("Spectrum of Gaussian (N=512)")
fig6[0].grid("true")
fig6[1].grid("true")
plt.show()
print(err)