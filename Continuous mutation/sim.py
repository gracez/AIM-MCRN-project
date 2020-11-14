import numpy as np
import matplotlib.pyplot as plt

###################################################
# simulates infection WITHOUT recurrent spillover 
# allows both gamma and beta to mutate            
#
# note: this implementation branches people in "generation order" instead of in "time order"
# (faster to run, and still gives the same result in the limit as outbreak_thresh -> infinity)
# for an implementation done in time order, see the function sim_class
#
# also note that if the infection goes extinct the time of extinction is well-defined, and either implementation computes it correctly.
###################################################
#
# input beta_0 is transmission rate of the initial case.
# input gamma_0 is recovery rate of the initial case. 
#
# optional input mu_1 represents mutation rate for beta. If unspecified, defaults to 0.
# optional input mu_2 represents mutation rate for gamma. If unspecified, defaults to 0.
#
###################################################
#
# returns 0, t_ext if extinction, where t_ext is time of extinction. 
# returns 1 if outbreak
#
###################################################

def sim(beta_0, gamma_0, mu_1=0, mu_2=0):
    
    outbreak_thresh = 100 #call it an outbreak if this many people are infected
    
    t = 0
    
    #store active cases as an array of triples
    #each triple looks like [beta, gamma, time of infection] and represents one infected person
    infecteds=[[beta_0, gamma_0, t]]
    
    N_infected = len(infecteds) #keep track of number of infected people
    
    t_ext = None #time of extinction, will fill this value and return it if the infection goes extinct
    
    while True: 
        
        #print("\n", "infecteds now = ", infecteds)
        
        #keeps track of indices of recovered people
        #will delete them from array of infecteds at the end of each while loop iteration
        recovereds = []
        
        for i, person in enumerate(infecteds): #go through and branch each person in array of infecteds

            #print("\n", "branching person", i)
            
            [beta, gamma, t] = person
            
            while True:
                                
                #interevent time
                #note this is interevent time for ONE person, NOT for entire population
                dt = np.random.exponential(scale = 1/(beta + gamma))
                new_t = t + dt #time of new event
                
                #pick which new event happens, transmission or recovery         
                ev = np.random.rand() #draw number from uniform distribution over [0, 1)
                prob_trans = beta / (beta+gamma)

                if (ev < prob_trans): #transmission
                    
                    #print("person", i, " is transmitting")
                    
                    #pick beta and gamma for new case
                    #mutation is a number drawn from normal distribution with std dev mu_1 or mu_2
                    #don't allow negative beta, gamma
                    mut1 = np.random.normal(loc=0.0, scale=mu_1)
                    new_beta = max(0, beta + mut1)
                    mut2 = np.random.normal(loc=0.0, scale=mu_2)              
                    new_gamma = max(0, gamma + mut2)
                    
                    #append new person to array of infecteds
                    infecteds.append([new_beta, new_gamma, new_t])
                    N_infected += 1
                    
                    if N_infected >= outbreak_thresh:
                        #print("\n", "outbreak!")
                        return 1

                    #print("new infection: ", [new_beta, new_gamma, new_t])

                else: #recovery
                    
                    #print("person", i, " recovers") 
                    
                    N_infected -= 1
 
                    if N_infected == 0:
                        #print("\n", "extinction")
                        return 0, new_t
                    
                    #mark index of this recovered person, in order to delete them later
                    #(we don't delete yet, because that will mess up the indexing of the for loop)
                    recovereds.append(i)
                    break #skip to next person

        #after each time we complete a round of branching everybody:
        #update the array of infecteds by deleting recoveries
        infecteds = np.delete(infecteds, recovereds, axis=0)