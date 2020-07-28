# AIM-MCRN-Summer-School-2020-Project-Zoonotic-Spillover
Zoonotic spillover research project for the COVID-19 summer school.

Mutation Subgroup.


TO DO: update this description!!



Simulation for branching process with random mutation. Each time the pathogen transmits to a new person, it mutates (both beta and gamma are allowed to mutate independently).

There are 2 main functions, "sim" and "recurrent_sim".

"sim" represents a process without recurrent spillover. It seeds the population with an initial case (presumably spilled over from a reservoir). It branches until either the infection goes extinct or surpasses a threshold. In the latter case, we call it an outbreak. 

"recurrent_sim" represents a process with random recurrent spillover at some rate alpha. It begins with a single founding spillover case and branches it. Meanwhile, more spillovers arrive, each giving birth to a new branching process tree. Overall, the infection either surpasses an outbreak threshold, or it goes extinct (which requires a moment in time when ALL spillover processes have gone extinct - i.e. there are 0 infected people at that time). 

To dos: 
    - Make more graphs/visualizations.
    - Play with mutation rates.
            - Play with unequal mutation rates for beta and gamma
    - Play with spillover rate.
    - Include backwards transmission.
    - Try with a finite population.
    - Do a sanity check for the implementation of interevent time against Linda's slides.
    - Do a sanity check for the bulk recurrent simulation against the Voinson 2018 paper.