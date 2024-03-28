# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math as m



def acceler(we,wi,x,y,z,t): #el temps hauria de ser l'interval pel qual volem simular
 # pero per ara el temps es un sol instant i ja ho arreglare dspres perq no se treballar amb python

    wx = we
    wy = wi*m.cos(we*t)
    wz = wi*m.sin(we*t)
    wv = np.array([wx, wy, wz])
    
    dwx = 0
    dwy = -wi*we*np.sin(we*t)
    dwz = wi*we*np.cos(we*t)
    dwv = np.array([dwx, dwy, dwz])
    
    rx = x*m.cos(wi*t)+z*m.sin(wi*t)
    ry = y*m.cos(we*t)+x*m.sin(we*t)*m.sin(wi*t)-z*m.sin(we*t)*m.cos(wi*t)
    rz = y*m.sin(we*t)-x*m.cos(we*t)*m.sin(wi*t)+z*m.cos(we*t)*m.cos(wi*t)
    rv = np.array([rx, ry, rz])
    
    wr = np.cross(wv,rv)
    ac =  np.cross(wv,wr)
    at =  np.cross(dwv,rv)
    
    
    atot = ac+at
    return ac, at, atot

###############################################################################
###############################################################################


# Define the parameters for our simulation:
we = 10*2*m.pi/60; #angular velocity of the exterior ring, in rad/s
wi = 5*2*m.pi/60; #angular velocity of the inner ring, in rad/s
x = 0.001; y = 0.001; z = 0.001; #position were we want to calculate the accelerations, in m
T = 120; #time of the simulations, in s

g = 9.8; #en valor absolut, al paper ho fan amb -9.8 pero crec q esta malament al paper

dt = 0.1; timevec = np.linspace(0,T,m.floor(T/dt)+1)
L = len(timevec)
Atotx = []; Aintx = []; Gintx = []; Aintx2 = []
Atoty = []; Ainty = []; Ginty = []; Ainty2 = []
Atotz = []; Aintz = []; Gintz = []; Aintz2 = []

taAx_int = []; taAy_int = []; taAz_int = []
taAx_int_nong = []; taAy_int_nong = []; taAz_int_nong = []

for i in range(L):
    t = timevec[i]
    ac, at, atot = acceler(we, wi, x, y, z, t)
    
    ax = atot[0]; ay = atot[1]; az = atot[2]
    
    # lab reference, no gravity
    Atotx = np.append(Atotx,ax)
    Atoty = np.append(Atoty,ay)
    Atotz = np.append(Atotz,az)
    
    axint = ax*m.cos(wi*t) + ay*m.sin(wi*t)*m.sin(we*t) - az*m.sin(wi*t)*m.cos(we*t)
    ayint = ay*m.cos(we*t) + az*m.sin(we*t)
    azint = ax*m.sin(wi*t) - ay*m.cos(wi*t)*m.sin(we*t) + az*m.cos(we*t)*m.cos(wi*t)
    
    # internal frame reference, no gravity
    Aintx = np.append(Aintx,axint)
    Ainty = np.append(Ainty,ayint)
    Aintz = np.append(Aintz,azint)
    
    tax = np.average(Aintx)
    tay = np.average(Ainty)
    taz = np.average(Aintz)
    
    taAx_int_nong = np.append(taAx_int_nong,tax)
    taAy_int_nong = np.append(taAy_int_nong,tay)
    taAz_int_nong = np.append(taAz_int_nong,taz)
    
    gxint = -g*m.sin(wi*t)*m.sin(we*t);     Gintx = np.append(Gintx,gxint)
    gyint = -g*m.cos(we*t);                 Ginty = np.append(Ginty,gyint)
    gzint = g*m.cos(wi*t)*m.sin(we*t);      Gintz = np.append(Gintz,gzint)
    
    # internal frame reference, with gravity
    Aintx2 = np.append(Aintx2,axint+gxint)
    Ainty2 = np.append(Ainty2,ayint+gyint)
    Aintz2 = np.append(Aintz2,azint+gzint)
    
    taxg = np.average(Aintx2)
    tayg = np.average(Ainty2)
    tazg = np.average(Aintz2)
    
    taAx_int = np.append(taAx_int,taxg)
    taAy_int = np.append(taAy_int,tayg)
    taAz_int = np.append(taAz_int,tazg)
    
    
fig,ax = plt.subplots() # lab reference accelerations, no gravity
ax.plot(timevec, Atotx, 'g')
ax.plot(timevec, Atoty, 'b')
ax.plot(timevec, Atotz, 'k')
plt.xlabel('time (s)')
plt.ylabel('acceleration ($m/s^2$)')
plt.title('Acceleration in laboratory frame of reference, w/o gravity')
plt.grid(True)
ax.legend('x''y''z')
plt.show()    


#Atoty2 = Atoty+g
#Atot = np.add(Atotx,Atoty2); Atot = np.add(Atot,Atotz)
#taA = np.average(Atot) #time averaged acceleration
#print(taA)


#fig2,ax2 = plt.subplots()
#ax2.plot(timevec, Atotx,'g')
#ax2.plot(timevec, Atoty2,'b')
#ax2.plot(timevec, Atotz,'k')
#ax2.plot(timevec, Atot,'r')
#plt.show()  



fig3,ax3 = plt.subplots() # internal accelerations, without gravity
ax3.plot(timevec, Aintx,'g')
ax3.plot(timevec, Ainty,'b')
ax3.plot(timevec, Aintz,'k')
#ax3.plot(timevec, Atot_int,'r')
plt.xlabel('time (s)')
plt.ylabel('acceleration ($m/s^2$)')
plt.title('Acceleration in sample frame of reference, w/o gravity')
ax3.legend('x''y''z')
plt.grid(True)
plt.show()
  
fig4,ax4 = plt.subplots() # internal accelerations, with gravity
ax4.plot(timevec, Aintx2,'g')
ax4.plot(timevec, Ainty2,'b')
ax4.plot(timevec, Aintz2,'k')
plt.xlabel('time (s)')
plt.ylabel('acceleration ($m/s^2$)')
plt.title('Acceleration in sample frame of reference, w/ gravity')
ax4.legend('x''y''z')
#ax4.plot(timevec, Atot_int,'r')
plt.grid(True)
plt.show()

taAx_int = np.array(taAx_int)
taAy_int = np.array(taAy_int)
taAz_int = np.array(taAz_int); 
taAtot_int = (taAx_int**2+taAy_int**2+taAz_int**2)**(.5)

fig5,ax5 = plt.subplots() # time averaged internal accelerations, with gravity
ax5.plot(timevec, taAx_int,'g')
ax5.plot(timevec, taAy_int,'b')
ax5.plot(timevec, taAz_int,'k')
ax5.plot(timevec, taAtot_int,'r')
plt.xlabel('time (s)')
plt.ylabel('acceleration ($m/s^2$)')
plt.title('Time averaged acceleration, w/ gravity')
ax5.legend(['x','y','z','total'])
plt.grid(True)
plt.show()

taAx_int_nong = np.array(taAx_int_nong)
taAy_int_nong = np.array(taAy_int_nong)
taAz_int_nong = np.array(taAz_int_nong); 
taAtot_int_nong = (taAx_int_nong**2+taAy_int_nong**2+taAz_int_nong**2)**(.5)

fig6,ax6 = plt.subplots() # time averaged internal accelerations, without gravity
ax6.plot(timevec, taAx_int_nong,'g')
ax6.plot(timevec, taAy_int_nong,'b')
ax6.plot(timevec, taAz_int_nong,'k')
ax6.plot(timevec, taAtot_int_nong,'r')
plt.xlabel('time (s)')
plt.ylabel('acceleration ($m/s^2$)')
plt.title('Time averaged acceleration, w/o gravity')
ax6.legend(['x','y','z','total'])
plt.grid(True)
plt.show()


max_effective_acc = max(taAtot_int[-500:-1])
print(max_effective_acc)
    
"""
Spyder Editor

This is a temporary script file.
"""

