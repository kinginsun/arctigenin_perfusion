;; 1. Based on: run101
;; 2. Description: 1-comp iv, linear elim
;; x1. Author: user
;; 3. Label:
$PROB BAYES TEST 

$DATA OLSID.CSV IGNORE=# 

$INPUT ID TIME OBS=DV

$PRED 
K = THETA(1) + ETA(1) 
F = 10*EXP(-K*TIME) 
IPRED = F 
Y = F + ERR(1) 

$THETA 10 ; K

$OMEGA 4 ; ETA(1)

$SIGMA 0.04 ; ERR(1)

$ESTIMATE MAXEVALS=0 POSTHOC

$TABLE ID TIME K IPRED ETA1 NOPRINT NOAPPEND ONEHEADER FILE=patab102.tab

