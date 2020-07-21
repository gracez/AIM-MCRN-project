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
## Set the Model Parameters and Initial Conditions
######################################################
# Rate of infection from infectious individuals
beta = 0.2
# Rate of infection from animal resevoir
tau = 0.0001
# Recovery rate
gamma = 0.1
# Waning rate of immunity
alpha = 0.0001

# Total population size
N = 200
# Initial infecteds
init = 5
# Initial conditions as a function of N and init
X0 = [N-init,init,0]

#%% ##################################################
## Set the Numerical Integration Parameters
######################################################
# Time step
dt = 0.01
# Starting time
tmin = 0
# Ending time
tmax = 100
# Integration mesh
t = np.arange(tmin,tmax,dt)
# Simulations
sim=3

Tend=np.zeros(len(t))
casetot=np.zeros(len(t))

for k in range(sim):
    
    i=np.zeros(len(t))
    r=np.zeros(len(t))
    s=np.zeros(len(t))
    tt=np.zeros(len(t))
    
    cases=init
    i[0]=init
    s[0]=N-init
    r[0]=0
    tt[0]=0
    j=0
    tot=N
    
    while (i[j]>=0 and j<=(len(t)-2)):
        omega=beta*s[j]*i[j]/tot+tau*s[j]+g*i[j]+alpha*r[j] #all transition rates
        
        u1=np.random.rand() #For interevent time
        u2=np.random.rand() #For the event that occurs
        
        tt[j+1]=tt[j]-np.log(u1)/omega #Interevent time
        
        event1=beta*i[j]*s[j]/tot/omega
        event2=event1+tau*s[j]/omega
        event3=event2+g*i[j]/omega
        event4=event3+alpha*r[j]/omega
        
        if (u2<=event1): #Transmission
            i[j+1]=i[j]+1
            s[j+1]=s[j]-1
            r[j+1]=r[j]
            cases=cases+1 #count total number of cases
            
        elif (u2>event1 and u2<=event2): #spillover
            i[j+1]=i[j]+1
            s[j+1]=s[j]-1
            r[j+1]=r[j]
            cases=cases+1 #count total number of cases
        
        elif (u2>event2 and u2<=event3): #recovery
            i[j+1]=i[j]-1
            r[j+1]=r[j]+1
            s[j+1]=s[j]
            
        elif (u2>event3 and u2<=event4): #waning immunity
            i[j+1]=i[j]
            r[j+1]=r[j]-1
            s[j+1]=s[j]+1
            
        else:
            i[j+1]=i[j]
            r[j+1]=r[j]
            s[j+1]=s[j]
        
        j=j+1
        tot=s[j]+i[j]+r[j]
    
    casetot[k]=cases #total number of cases for each simulation
    Tend[k]=tt[j] #time epidemic ends for each simulation
    
    if k==0:
        plt.step(tt[0:j],i[0:j],'r-',label='k1')
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')
        plt.grid(alpha=0.33)
        plt.xlim(0,200)
        
    if k==1:
        plt.step(tt[0:j],i[0:j],'b-',label='k2')
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')
        plt.grid(alpha=0.33)

        
print("Case1={},Case2={}".format(casetot[0],casetot[1]))
print()
print("Tend1={},Tend2={}".format(Tend[0],Tend[1]))

