#!/usr/bin/env python
# coding: utf-8

# In[415]:


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
beta = 0.09
# Rate of infection from animal resevoir
tau = .001
# Recovery rate
gamma = 0.1
# Waning rate of immunity
alpha = 0.01

B0=tau
B1=.9

# Total population size
N = 100
# Initial infecteds
init = 1
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
tmax = 365*4
# Integration mesh
t = np.arange(tmin,tmax,dt)

plt.rcParams["figure.figsize"] = (12,6)


# In[416]:



#%% ##################################################
## Define the Zoonotic ODE Model
######################################################
# A three-dimensional system (S,I,R) where:
# X[0] = dS
# X[1] = dI
# X[2] = dR
def dX_dt(X,t):
    
    return [-beta*X[0]*X[1]/N - tau*X[0] + alpha*X[2],
           beta*X[0]*X[1]/N + tau*X[0] - gamma*X[1],
           gamma*X[1] - alpha*X[2]]

#%% ##################################################
## Integrate the System and Store Output
######################################################
# Solve the system
ansode = odeint(dX_dt,X0,t)

# Save state variable dyanimcs
susode = ansode[:,0]
infode = ansode[:,1]
recode = ansode[:,2]

#%% ##################################################
## Plot Findings
######################################################
plt.plot(t,susode,label='susceptible',linewidth=2)
plt.plot(t,infode,label='infected',linewidth=2,linestyle='--')
plt.plot(t,recode,label='recovered',linewidth=2,linestyle=':')
plt.xlabel('time',fontsize=13)
plt.ylabel('population',fontsize=13)
plt.legend()
plt.grid(alpha=0.33)
plt.show()


# In[417]:



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

#%% ##################################################
## Integrate the System and Store Output
######################################################
# Solve the system
anssquare = odeint(dXsquare_dt,X0,t)
anscos = odeint(dXcos_dt,X0,t)

# Save state variable dyanimcs
sussq = anssquare[:,0]
infsq = anssquare[:,1]
recsq = anssquare[:,2]

suscos = anscos[:,0]
infcos = anscos[:,1]
reccos = anscos[:,2]

#%% ##################################################
## Plot Findings
######################################################
plt.plot(t,sussq,label='susceptible',linewidth=2)
plt.plot(t,infsq,label='infected',linewidth=2,linestyle='--')
plt.plot(t,recsq,label='recovered',linewidth=2,linestyle=':')
plt.xlabel('time',fontsize=13)
plt.ylabel('population',fontsize=13)
plt.legend()
plt.grid(alpha=0.33)
plt.show()

plt.plot(t,suscos,label='susceptible',linewidth=2)
plt.plot(t,infcos,label='infected',linewidth=2,linestyle='--')
plt.plot(t,reccos,label='recovered',linewidth=2,linestyle=':')
plt.xlabel('time',fontsize=13)
plt.ylabel('population',fontsize=13)
plt.legend('best')
plt.grid(alpha=0.33)
plt.show()

plt.plot(t,B0*(1+-B1*np.cos((2*np.pi/365)*t)))


# In[371]:



#CTMC
#%% ##################################################
## Set the Numerical Integration Parameters
######################################################
# Simulations
sim=3

Tend=np.zeros(len(t))
casetot=np.zeros(len(t))

plt.plot(t,infode)

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
        omega=beta*s[j]*i[j]/tot+tau*s[j]+gamma*i[j]+alpha*r[j] #all transition rates
        
        u1=np.random.rand() #For interevent time
        u2=np.random.rand() #For the event that occurs
        
        tt[j+1]=tt[j]-np.log(u1)/omega #Interevent time
        
        event1=beta*i[j]*s[j]/tot/omega
        event2=event1+tau*s[j]/omega
        event3=event2+gamma*i[j]/omega
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
        plt.xlim(0,400)
        
    if k==1:
        plt.step(tt[0:j],i[0:j],'b-',label='k1')
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')
        plt.grid(alpha=0.33)

        
print("Case1={},Case2={}".format(casetot[0],casetot[1]))
print()
print("Tend1={},Tend2={}".format(Tend[0],Tend[1]))


# In[372]:



#parameters
Tend=np.zeros(len(t)) #initializing a vector of zeros
casetot=np.zeros(len(t))
count=0
sim=3

for k in range(sim):
    
    i=np.zeros(len(t))
    r=np.zeros(len(t))
    s=np.zeros(len(t))
    
    cas=init
    i[0]=init
    s[0]=N-init
    r[0]=0
    tot=N
    j=0
    tim=0
    
    while (i[j]>=0 and j<len(t)-1):
        u=np.random.rand() #uniform random number
        
        #tau=B0*(1+B1*np.cos(2*np.pi/365*tim*dt))
        
        if np.cos(2*np.pi/365*tim*dt)<0:
            tau=B0*(1+B1*-1)
            
        elif np.cos(2*np.pi/365*tim*dt)>0:
            tau=B0*(1+B1*1) 
        else:
            tau=B0*(1)
           
        
        ev1=beta*i[j]*s[j]/tot*dt
        ev2=ev1+tau*s[j]*dt 
        ev3=ev2+gamma*i[j]*dt
        ev4=ev3+alpha*r[j]*dt
        
        if ev2>1:
            count=count+1 #should always be zero
        
        if (u<=ev1): #transmission
            i[j+1]=i[j]+1
            s[j+1]=s[j]-1
            r[j+1]=r[j]
            cas=cas+1 #count total number of cases
        
        elif (u>ev1 and u<=ev2): #spillover
            i[j+1]=i[j]+1
            r[j+1]=r[j]
            s[j+1]=s[j]-1
            
        elif (u>ev2 and u<=ev3): #recovery
            i[j+1]=i[j]-1
            r[j+1]=r[j]+1
            s[j+1]=s[j]
            
        elif (u>ev3 and u<=ev4): #waned immunity
            i[j+1]=i[j]
            r[j+1]=r[j]-1
            s[j+1]=s[j]+1
        
        else: #no change
            i[j+1]=i[j]
            r[j+1]=r[j]
            s[j+1]=s[j]
             
        j=j+1
        tim=tim+1;
        tot=s[j]+i[j]+r[j]
        t[j]=t[j-1]+dt
    
    casetot[k]=cas #total cases for each simulation
    Tend[k]=t[j] #time epidemic ends for each simulation
    
    if k==0:
        plt.step(t,i,'m-.',label='k1')
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')
      
    elif k==1:
        plt.step(t,i,'g-.',label='k2') 
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')
        
    elif k==2:
        plt.step(t,i,'b-.',label='k3')
        plt.xlabel('Day')
        plt.ylabel('Infectious Individuals')

plt.plot(t,infsq,'k-')
plt.title("DTMC SIR")
plt.legend(loc='best')
plt.grid(alpha=0.33)
plt.show()



print("Case1={}, Case2={}, Case3={}".format(casetot[0], casetot[1], casetot[2]))
print()
print("Prob>1 ={}".format(count))


# In[ ]:





# In[418]:


N = 100
#Euler Maruyama Stochastic Solution
S=np.zeros(len(t))
I=np.zeros(len(t))
R=np.zeros(len(t))
Tau=np.zeros(len(t))
m=3

for i in range(m):
    print(i)
    S[0]=N-init
    I[0]=init
    R[0]=0
    Tau[0]=0
    j=0
    
    
    if (N<101 and N>99):
        while (j<=(len(t)-3)):
            if (I[j]>=0):
                #f=B0*(1+B1*np.cos((2*np.pi/365)*t[j]))+Tau[j]
                
                if np.cos(2*np.pi/365*t[j])<0:
                    f=B0*(1+B1*-1)+Tau[j]

                elif np.cos(2*np.pi/365*t[j])>0:
                    f=B0*(1+B1*1) +Tau[j]
                else:
                    f=B0*(1)+Tau[j]

                S[j+1]=S[j]+(-(beta*S[j]*I[j]/N)-(f*S[j])+alpha*R[j])*dt
                I[j+1]=I[j]+(beta*S[j]*I[j]/N+(f*S[j])-gamma*I[j])*dt
                R[j+1]=R[j]+(gamma*I[j]-alpha*R[j])*dt
                Tau[j+1]=Tau[j]+0.0001*np.random.randn()*np.sqrt(dt)
                j=j+1
                N=S[j]+I[j]+R[j]

            else:
                I[j]=0

                if np.cos(2*np.pi/365*t[j])<0:
                    f=B0*(1+B1*-1)+Tau[j]

                elif np.cos(2*np.pi/365*t[j])>0:
                    f=B0*(1+B1*1) +Tau[j]
                else:
                    f=B0*(1)+Tau[j]

                S[j+1]=S[j]+(-(beta*S[j]*I[j]/N)-(f*S[j])+alpha*R[j])*dt
                I[j+1]=I[j]+(beta*S[j]*I[j]/N+(f*S[j])-gamma*I[j])*dt
                R[j+1]=R[j]+(gamma*I[j]-alpha*R[j])*dt
                Tau[j+1]=Tau[j]+0.0001*np.random.randn()*np.sqrt(dt)
                j=j+1
                N=S[j]+I[j]+R[j]
        
    else:
        pass
        
    plt.plot(t[0:j],I[0:j])
    plt.plot(t,infsq,'k--',label=('det'))
    plt.xlabel('Days')
    plt.ylabel('Infectious')


# In[ ]:




