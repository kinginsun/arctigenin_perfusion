;; 1. Based on: run107
;; 2. Description: single function
;; x1. Author: user
;; 3. Label:
$PROB BAYES TEST 

$DATA data_100.csv IGNORE=# 

$INPUT ID TIME OBS=DV

$PRED 
K = THETA(1) + ETA(1) 
F = 1000*EXP(-K*TIME) 
IPRED = F 
Y = F + ERR(1) 

$THETA 1 ; K

$OMEGA 0.1 ; ETA(1)

$SIGMA 0.04 ; ERR(1)

$ESTIMATE MAXEVALS=9999 METHOD=0

$TABLE ID TIME K IPRED NOPRINT NOAPPEND ONEHEADER FILE=patab108.tab

