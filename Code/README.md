## Overview
`util.R` contains the following functions
* LDM_test:
* LINDA_test
* LinDA_test_unadj
* MaAsLin2_test
* MaAsLin2_test_unadj 
* NBZIMM_test
* NBZIMM_test_unadj
* NBZIMM_indiv_test
* ZIBR_test
* QRank_test
* minP_combine_test
* cauchy_combine_test
* zero_inf_ind_rank_score_test: this is the main function of ZINQ-L, it will call `minP_combine_test` and `cauchy_combine_test` to aggregate tests from multiple quantile levels.

Note: 
* `XXX_test_unadj` is the unadjusted version of `XXX_test`, which only includes the intercept term without adjusting for other covariates. For `XXX_test`, the list of adjusting covariates is extracted from our example data. If you want to run on your own data, you need to modify the covariates name in the function. 

* A bunch of messages will pop up when running some existing methods, e.g., MaAsLin2 will create a temporary folder inside the repo with names `tmp_{job.id}` where `job.id` is the user specified id that controls the sim setting (default: 1).
