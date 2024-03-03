## Overview

* `sim2_complete.R`: simulation scenario 2 - adjusted analysis on OTU table under complete null or alternative. The output will be a list of vectors representing power of tests, where each element of the list corresponds to one simulation run, and each vector contains power of the tests listed in the following section. 

* `sim3_partial.R`: simulation scenario 3 - adjusted analysis on OTU table with partial null and alternative. The output will be a list of data frames. Each data frame contains the FDP (false disvcovery proportion) and TPR (True positive rate) of all the tests. 


## Tests
* lm
* lmm
* ZINQ-L
  * independent correlation structure (ZINQ-L reduce to ZINQ) across 5 levels
  * dependent correlation structure (ZINQ-L) across 5 levels
  * minP_combine_test that use 5 quantile levels or clip to 3 quantile levels (0.25, 0.5, 0.75)
  * cauchy_combine_test with equal/unequal weights
* LDM (comment out to reduce computation burden)
* LINDA
* MaAsLin2
* under NBZIMM umbrella
  * NB
  * ZINB
  * ZIG_count
  * ZIG_prop


Note: 
* Please refer to `example.R` for simulation scenario 1 - unadjusted analysis on individual taxa.  

* A bunch of messages will pop up when running some existing methods. Users can create a folder 'Scratch' to save intermediate output, e.g., MaAsLin2 will create a temporary folder inside the 'Scratch' folder with names `tmp_{job.id}` where `job.id` is the user specified id that controls the sim setting (default: 1). To suppress all the messages, use the `hush` function defined in util.R. 
