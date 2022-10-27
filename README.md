# BDS_deterministic_map

- main script is called `ODE_SCM.R`



The last concern was more crucial, namely if we should multiply D(t), the backwards probabilities/conditional likelihoods, into the forwards ODEs. Just to repeat and clarify the idea, we want to sample the current rate category at time t+\Delta \delta with \Delta t < \Delta \delta (we have to sample at pre-defined times along the branches every \Delta\delta time steps). We know the forward probability F(t) and the backwards probability D(t+\Delta \delta). The big question is how do we get from F(t) to F(t+\Delta \delta). To find the probability for sampling the rate categories, we have 3 different solutions, which I denoted A, B and C. (A) we set up F(t), use the ODEs that multiply with D(t), and integrate them until t+\Delta \delta; (B) we initialize F(t), then integrate it forward in time until t+\Delta \delta, and then multiply it with D(t+\Delta \delta); or (C) we randomly pick the category at time t with probabilities F(t), then initialize F_i(t)=1 and F_j(t)=0 for all j!=i, next we integrate F(t) forward in time, and finally multiply it with D(t+\Delta \delta). I have implemented these three different ways into the attached R scripts. I have also set up scripts to estimate the ancestral state probabilities using diversitree, our RevBayes ancestral state monitor, and our RevBayes stochastic character mapping algorithm.

The results of these test show a few things. First, our RevBayes implementation gives the same answer as diversitree. I think this really validates that our implementation is correct (differences are due to sampling uncertainty). Second, only algorithms B and C are correct, but not algorithm A.

 Equation 5 and 6 should be implemented
