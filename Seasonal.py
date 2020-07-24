#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#%% ##################################################
## Import Libraries
######################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#%% ##################################################
## Define the Zoonotic ODE Model
######################################################
# A three-dimensional system (S,I,R) where:
# X[0] = dS
# X[1] = dI
# X[2] = dR
def dXsquare_dt(X,t):
    
    if np.cos(2*np.pi/365*t)<0:
        tau=B0*(1+B1*-1)  
    elif np.cos(2*np.pi/365*t)>0:
        tau=B0*(1+B1*1) 
    else:
        tau=B0*(1)
    
    return [-beta*X[0]*X[1]/N - tau*X[0] + alpha*X[2],
           beta*X[0]*X[1]/N + tau*X[0] - gamma*X[1],
           gamma*X[1] - alpha*X[2]]

def dXcos_dt(X,t):
    
    tau=B0*(1+B1*np.cos((2*np.pi/365)*t))
    
    return [-beta*X[0]*X[1]/N - tau*X[0] + alpha*X[2],
           beta*X[0]*X[1]/N + tau*X[0] - gamma*X[1],
           gamma*X[1] - alpha*X[2]]

def dXseas_dt(X,t):
    
    tau=B0/(1+np.exp(B1*(t-79)**2))
    
    return [-beta*X[0]*X[1]/N - tau*X[0] + alpha*X[2],
           beta*X[0]*X[1]/N + tau*X[0] - gamma*X[1],
           gamma*X[1] - alpha*X[2]]

#%% ##################################################
## Set the Model Parameters and Initial Conditions
######################################################
# Rate of infection from infectious individuals
beta = 0.57
# Rate of infection from animal resevoir
tau2 = 0.001
# Recovery rate
gamma = 0.034
# Waning rate of immunity
alpha = 0.01

B0=1
B1=0.001

# Total population size
N = 100
# Initial infecteds
init = 0
# Initial conditions as a function of N and init
X0 = [N-init,init,0]

#%% ##################################################
## Set the Numerical Integration Parameters
######################################################
# Time step
dt = 0.005
# Starting time
tmin = 0
# Ending time
tmax = 365
# Integration mesh
t = np.arange(tmin,tmax,dt)

#%% ##################################################
## Integrate the System and Store Output
######################################################
# Solve the system
ans1 = odeint(dXsquare_dt,X0,t)
ans2 = odeint(dXcos_dt,X0,t)
# Weird-shape wave gets special treatment!
addyrs = 2
ans3 = odeint(dXseas_dt,X0,t)
ts = t

for i in range(addyrs):
    X0 = ans3[-1]
    ans3 = np.concatenate((np.delete(ans3,np.array([len(ans3)-1]),axis=0),odeint(dXseas_dt,X0,t)),axis=0)
    ts = np.concatenate((np.delete(ts,np.array([len(ts)-1])),t+ts[-1]),axis=0)
    

# Save state variable dyanimcs
sus1 = ans1[:,0]
inf1 = ans1[:,1]
rec1 = ans1[:,2]

sus2 = ans2[:,0]
inf2 = ans2[:,1]
rec2 = ans2[:,2]

sus3 = ans3[:,0]
inf3 = ans3[:,1]
rec3 = ans3[:,2]

#%% ##################################################
## Plot Findings
######################################################
# plt.plot(t,sus1,label='susceptible',linewidth=2)
# plt.plot(t,inf1,label='infected',linewidth=2,linestyle='--')
# plt.plot(t,rec1,label='recovered',linewidth=2,linestyle=':')
# plt.xlabel('days',fontsize=13)
# plt.ylabel('population',fontsize=13)
# plt.legend()
# plt.grid(alpha=0.33)
# plt.show()

# plt.plot(t,sus2,label='susceptible',linewidth=2)
# plt.plot(t,inf2,label='infected',linewidth=2,linestyle='--')
# plt.plot(t,rec2,label='recovered',linewidth=2,linestyle=':')
# plt.xlabel('days',fontsize=13)
# plt.ylabel('population',fontsize=13)
# plt.legend()
# plt.grid(alpha=0.33)
# plt.show()

plt.plot(ts,sus3,label='susceptible',linewidth=3)
plt.plot(ts,inf3,label='infected',linewidth=3,linestyle='--')
plt.plot(ts,rec3,label='recovered',linewidth=3,linestyle=':')
plt.xlabel('days',fontsize=14)
plt.ylabel('population',fontsize=14)
plt.title('Summer spillovers',fontsize=16)
plt.legend(fontsize=12)
plt.grid(alpha=0.33)
plt.show()

